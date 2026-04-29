#!/usr/bin/env python3
"""For each `no_path` failure case, check whether Neo4j has a path in the
original pathway. Distinguishes real generation bugs from genuine v96
pathway changes.

Methodology per failure case:
  1. Load the failing pair (gene, keyoutput, pathway) from the failures TSV.
  2. Build the pathway's reaction-flow graph from Neo4j: edges are
     (entity → reaction) for inputs/catalysts, (reaction → entity) for
     outputs. Catalysts and regulators are "consumes" relationships in this
     model — perturbing the catalyst should affect the reaction's output.
  3. Map the gene's stable IDs (any modification/compartment form) and the
     key-output's stable ID into nodes in that graph.
  4. BFS from gene-entities to keyoutput-entity.
  5. Categorize:
       - bug_candidate: Neo4j has a path; the logic network missed it.
       - pathway_changed: Neo4j has no path either; the curator's
         expectation was based on v86 connectivity that v96 has dropped.

Output: a categorized TSV alongside a summary count.
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.argument_parser import logger  # noqa: E402
from src.neo4j_connector import get_graph  # noqa: E402

MP_BIOPATH_DIR = Path("/home/awright/gitroot/mp-biopath-pathways")


def load_pathway_lookup() -> dict[str, str]:
    df = pd.read_csv("/tmp/mpbio_pathways.tsv", sep="\t")
    return dict(zip(df["pathway_name"], df["id"]))  # name → R-HSA-id


def build_neo4j_pathway_graph(graph, pathway_id: str) -> dict[str, list[str]]:
    """Forward adjacency: entity_stid → list of reachable_entity_stids by stepping through one reaction.

    For each reaction in the pathway, every input/catalyst/regulator entity
    is an edge source pointing at every output entity. Models 'perturbing
    this entity affects the production of these outputs' — the same flow
    the logic network represents.
    """
    rows = graph.run(
        """
        MATCH (p:Pathway {stId: $pid})-[:hasEvent*]->(r:ReactionLikeEvent)
        OPTIONAL MATCH (r)-[:input]->(in_pe:PhysicalEntity)
        WITH r, collect(DISTINCT in_pe.stId) AS inputs
        OPTIONAL MATCH (r)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(cat_pe:PhysicalEntity)
        WITH r, inputs, collect(DISTINCT cat_pe.stId) AS catalysts
        OPTIONAL MATCH (r)-[:regulatedBy]->(:Regulation)-[:regulator]->(reg_pe:PhysicalEntity)
        WITH r, inputs, catalysts, collect(DISTINCT reg_pe.stId) AS regulators
        OPTIONAL MATCH (r)-[:output]->(out_pe:PhysicalEntity)
        WITH r, inputs, catalysts, regulators, collect(DISTINCT out_pe.stId) AS outputs
        RETURN r.stId AS rxn, inputs, catalysts, regulators, outputs
        """,
        pid=pathway_id,
    ).data()

    adj: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        producers = set(row["inputs"]) | set(row["catalysts"]) | set(row["regulators"])
        producers.discard(None)
        outputs = [o for o in row["outputs"] if o]
        for src in producers:
            adj[src].extend(outputs)
    return dict(adj)


def expand_to_decomposition_leaves(graph, stids: list[str]) -> set[str]:
    """For each input stId, also include every entity it decomposes into via
    hasComponent / hasMember / hasCandidate (up to 5 levels deep). Lets a
    perturbation on a Complex-form be matched against entities in any
    sub-decomposition the logic network might have used.
    """
    if not stids:
        return set()
    rows = graph.run(
        """
        UNWIND $stids AS sid
        MATCH (root:PhysicalEntity {stId: sid})
        OPTIONAL MATCH (root)-[:hasComponent|hasMember|hasCandidate*0..5]->(child:PhysicalEntity)
        RETURN sid, collect(DISTINCT child.stId) AS children
        """,
        stids=stids,
    ).data()
    out: set[str] = set(stids)
    for r in rows:
        for c in r["children"]:
            if c:
                out.add(c)
    return out


def gene_to_stids(graph, gene_name: str) -> list[str]:
    rows = graph.run(
        """
        MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE $gene IN re.geneName
        RETURN COLLECT(DISTINCT pe.stId) AS stids
        """,
        gene=gene_name,
    ).data()
    return rows[0]["stids"] if rows else []


def reachable_in_graph(adj: dict[str, list[str]], sources: set[str]) -> set[str]:
    visited = set(sources)
    frontier = list(sources)
    while frontier:
        nxt = []
        for u in frontier:
            for v in adj.get(u, ()):
                if v not in visited:
                    visited.add(v)
                    nxt.append(v)
        frontier = nxt
    return visited


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--failures-tsv",
        default="validation_results/2026-04-29_mpbio_per_pathway_failures.tsv",
        help="Failures TSV from validate-against-mpbiopath.py",
    )
    ap.add_argument(
        "--out",
        default="validation_results/2026-04-29_no_path_neo4j_check.tsv",
        help="Output: per-case TSV with bug_candidate / pathway_changed / both / no_path classification",
    )
    args = ap.parse_args()

    failures = pd.read_csv(args.failures_tsv, sep="\t")
    no_path = failures[failures["category"] == "no_path"].copy()
    print(f"Total no_path cases: {len(no_path)}")

    pathway_lookup = load_pathway_lookup()
    graph = get_graph()

    # Cache per pathway: adjacency graph
    adj_cache: dict[str, dict[str, list[str]]] = {}
    # Cache per gene: decomposed stids
    gene_cache: dict[str, set[str]] = {}

    results = []
    bug_candidate = 0
    pathway_changed = 0
    cant_resolve = 0

    for _, row in no_path.iterrows():
        pname = row["pathway"]
        gene = row["gene"]
        ko = str(row["key_output"])

        pid = pathway_lookup.get(pname)
        if not pid:
            results.append({**row.to_dict(), "neo4j_status": "no_pathway_id"})
            cant_resolve += 1
            continue

        if pid not in adj_cache:
            adj_cache[pid] = build_neo4j_pathway_graph(graph, pid)
        adj = adj_cache[pid]

        if gene not in gene_cache:
            stids = gene_to_stids(graph, gene)
            decomposed = expand_to_decomposition_leaves(graph, stids)
            gene_cache[gene] = decomposed
        gene_stids = gene_cache[gene]

        target_stid = f"R-HSA-{ko}"

        if not gene_stids:
            results.append({**row.to_dict(), "neo4j_status": "gene_not_in_neo4j"})
            cant_resolve += 1
            continue

        reachable = reachable_in_graph(adj, gene_stids)
        if target_stid in reachable:
            status = "bug_candidate"
            bug_candidate += 1
        else:
            status = "pathway_changed"
            pathway_changed += 1
        results.append({**row.to_dict(), "neo4j_status": status})

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, sep="\t", index=False)
    print()
    print("=== Summary ===")
    print(f"  bug_candidate:    {bug_candidate} (Neo4j has a path; logic network missed it)")
    print(f"  pathway_changed:  {pathway_changed} (Neo4j has no path either; v86→v96 change)")
    print(f"  could not resolve: {cant_resolve}")
    if bug_candidate + pathway_changed > 0:
        rate = bug_candidate / (bug_candidate + pathway_changed) * 100
        print(f"  bug rate: {rate:.1f}% of resolvable no_path cases")
    print(f"Saved per-case classification to {args.out}")


if __name__ == "__main__":
    main()
