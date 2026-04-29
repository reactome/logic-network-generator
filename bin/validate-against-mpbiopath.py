#!/usr/bin/env python3
"""Validate generated logic networks against the MP-BioPath curator dataset.

For each pathway in the MP-BioPath validation set:
1. Load the regenerated logic_network.csv and stid_to_uuid_mapping.csv.
2. For each (perturbation, key-output) pair from the curator predictions,
   apply Boolean propagation and compare the predicted state to the curator's.
3. Aggregate accuracy and a confusion matrix per pathway and overall.

Propagation lattice: {0=down, 1=normal, 2=up}.

Per-node update:
- VR node (a UUID that appears in reaction_id_map.csv): the activity is the
  MIN of (a) each input/catalyst/positive-regulator incoming source state
  and (b) the INVERTED state of each negative-regulator incoming source.
- Entity node, normal case: the MAX of producer-VR activities (OR over
  producers, since output edges have and_or='or' when multiple producers).
- Entity node receiving 'assembly' edges (root-input complex): MIN of leaf
  states (assembly is AND).
- Entity node receiving a 'dissociation' edge (terminal-output leaf): the
  source complex's state.
- Pinned (perturbed) nodes hold their pinned state across iterations.

Loops: synchronous update with a 50-iteration cap. With three discrete
states the system either reaches a fixed point or starts oscillating; we
return the final synchronous state and accept that oscillating cases are
inherently uncertain.
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
DOWN, NORMAL, UP = 0, 1, 2
MAX_ITERATIONS = 50


def invert(state: int) -> int:
    return {DOWN: UP, NORMAL: NORMAL, UP: DOWN}[state]


def find_pathway_output_dir(output_root: Path, pathway_id: str) -> Path | None:
    """Locate the regenerated pathway dir given a numeric or R-HSA pathway id."""
    pid = pathway_id.replace("R-HSA-", "")
    for d in sorted(output_root.iterdir()):
        if d.is_dir() and d.name.endswith(f"R-HSA-{pid}"):
            return d
    return None


def load_network(pathway_dir: Path) -> pd.DataFrame:
    return pd.read_csv(pathway_dir / "logic_network.csv")


def build_stid_to_uuids(pathway_dir: Path) -> dict[str, list[str]]:
    """stable_id → list of UUIDs (multiple if entity at multiple positions)."""
    mapping = pd.read_csv(pathway_dir / "stid_to_uuid_mapping.csv")
    out: dict[str, list[str]] = defaultdict(list)
    for _, row in mapping.iterrows():
        out[str(row["stable_id"])].append(str(row["uuid"]))
    return out


def gene_name_to_stids(graph, gene_names: list[str]) -> dict[str, list[str]]:
    """Map gene names to all PhysicalEntity stable IDs whose reference entity has that gene name."""
    rows = graph.run(
        """
        UNWIND $names AS gene
        MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE gene IN re.geneName
        RETURN gene AS gene, COLLECT(DISTINCT pe.stId) AS stids
        """,
        names=gene_names,
    ).data()
    return {r["gene"]: r["stids"] for r in rows}


def build_incoming_index(network: pd.DataFrame) -> tuple[dict[str, list[tuple]], set[str]]:
    """Per-target list of (source, edge_type, pos_neg) tuples + universe of UUIDs.

    Built once per pathway via numpy iteration; reused across every
    perturbation's propagation run.
    """
    incoming: dict[str, list[tuple]] = defaultdict(list)
    all_uuids: set[str] = set()
    sources = network["source_id"].values
    targets = network["target_id"].values
    edge_types = network["edge_type"].values
    pos_negs = network["pos_neg"].values
    for s, t, et, pn in zip(sources, targets, edge_types, pos_negs):
        all_uuids.add(s)
        all_uuids.add(t)
        incoming[t].append((s, et, pn))
    return incoming, all_uuids


def propagate(
    incoming: dict[str, list[tuple]],
    all_uuids: set[str],
    pinned: dict[str, int],
) -> dict[str, int]:
    """Run Boolean propagation until fixed point or MAX_ITERATIONS."""
    state = {u: pinned.get(u, NORMAL) for u in all_uuids}

    for _ in range(MAX_ITERATIONS):
        new_state = dict(state)
        for uuid_, edges in incoming.items():
            if uuid_ in pinned:
                continue  # perturbed nodes stay pinned
            and_contribs: list[int] = []
            or_contribs: list[int] = []
            for source, etype, pn in edges:
                src_state = state.get(source, NORMAL)
                if etype == "regulator" and pn == "neg":
                    and_contribs.append(invert(src_state))
                elif etype in {"input", "catalyst", "regulator", "assembly"}:
                    and_contribs.append(src_state)
                elif etype in {"output", "dissociation"}:
                    or_contribs.append(src_state)

            if and_contribs and or_contribs:
                new_state[uuid_] = max(min(and_contribs), max(or_contribs))
            elif and_contribs:
                new_state[uuid_] = min(and_contribs)
            elif or_contribs:
                new_state[uuid_] = max(or_contribs)
        if new_state == state:
            break
        state = new_state
    return state


def predict(incoming, all_uuids, gene_to_uuids, key_output_uuids, perturbation_gene, direction):
    """Run one perturbation; return predicted states for each key-output uuid."""
    pinned: dict[str, int] = {}
    for u in gene_to_uuids.get(perturbation_gene, []):
        pinned[u] = direction
    final = propagate(incoming, all_uuids, pinned)
    return {ko: max((final.get(u, NORMAL) for u in uuids), default=NORMAL)
            for ko, uuids in key_output_uuids.items()}


def parse_perturbation_columns(curator_df: pd.DataFrame) -> list[tuple[str, int]]:
    """Columns like 'TP53_0' / 'TP53_2' → list of (gene, direction)."""
    out = []
    for col in curator_df.columns:
        if col in {"key_output", "control"}:
            continue
        if "_" in col:
            gene, suffix = col.rsplit("_", 1)
            if suffix in {"0", "2"}:
                out.append((gene, int(suffix), col))
    return out


_KEY_OUTPUT_ALIASES = ("key_output", "key output", "key_outout", "key outout")


def build_adjacency(network: pd.DataFrame) -> dict[str, list[str]]:
    """Forward adjacency list (source → list of targets) using numpy arrays.

    iterrows() is unusable on million-row DataFrames; zipping the underlying
    numpy arrays is two orders of magnitude faster.
    """
    adj: dict[str, list[str]] = defaultdict(list)
    for s, t in zip(network["source_id"].values, network["target_id"].values):
        adj[s].append(t)
    return adj


def reachable_from(adj: dict[str, list[str]], sources: set[str]) -> set[str]:
    """BFS over a precomputed adjacency dict."""
    if not sources:
        return set()
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


def categorize_failure(
    predicted: int,
    expected: int,
    perturbed_uuids: list[str],
    keyoutput_uuids: list[str],
    reachable: set[str],
) -> str:
    """Classify a failing test case so we can tell network bugs from propagator limits.

    - pass: predicted == expected (not a failure)
    - gene_not_in_network: the perturbed gene has no UUIDs in this network
    - keyoutput_not_in_network: the key-output entity has no UUIDs in this network
    - no_path: both endpoints exist but no directed path connects them
    - false_positive_change: predicted a perturbation but curator expected normal
    - propagator_missed: path exists, perturbed end was perturbed, but propagation
      didn't carry the change to the key output (or carried the wrong direction)
    """
    if predicted == expected:
        return "pass"
    if not perturbed_uuids:
        return "gene_not_in_network"
    if not keyoutput_uuids:
        return "keyoutput_not_in_network"
    if not any(ko in reachable for ko in keyoutput_uuids):
        return "no_path"
    if expected == 1 and predicted != 1:
        return "false_positive_change"
    return "propagator_missed"


def validate_one_pathway(pathway_dir: Path, pathway_name: str, pathway_dbid: str, graph) -> dict:
    """Validate one pathway, return per-pathway metrics."""
    curator_path = MP_BIOPATH_DIR / "reactome_curator_predictions" / f"{pathway_name}_reactome_curator_results.tsv"
    if not curator_path.exists():
        return {"status": "no_curator_file", "name": pathway_name}

    # Lenient parser tolerates trailing tabs and other small formatting drift
    # (some files have stray empty cells on certain rows).
    curator = pd.read_csv(curator_path, sep="\t", engine="python", on_bad_lines="warn")
    # Normalize the key-output column name to a single canonical "key_output"
    for alias in _KEY_OUTPUT_ALIASES:
        if alias in curator.columns:
            if alias != "key_output":
                curator = curator.rename(columns={alias: "key_output"})
            break
    else:
        return {"status": f"no key-output column (saw {list(curator.columns)[:3]}...)", "name": pathway_name}
    network = load_network(pathway_dir)
    stid_to_uuids = build_stid_to_uuids(pathway_dir)

    # Resolve key outputs (numeric dbId → list of UUIDs)
    key_output_uuids: dict[str, list[str]] = {}
    for ko in curator["key_output"].astype(str):
        key_output_uuids[ko] = stid_to_uuids.get(f"R-HSA-{ko}", [])

    # Resolve perturbation genes
    perturbations = parse_perturbation_columns(curator)
    gene_names = sorted({g for g, _, _ in perturbations})
    gene_to_stids = gene_name_to_stids(graph, gene_names)
    gene_to_uuids: dict[str, list[str]] = {}
    for g, stids in gene_to_stids.items():
        uuids = []
        for sid in stids:
            uuids.extend(stid_to_uuids.get(sid, []))
        gene_to_uuids[g] = uuids

    # Build adjacency and incoming-edge index once per pathway and reuse
    # them across every perturbation's reachability check and propagation.
    adj = build_adjacency(network)
    incoming, all_uuids = build_incoming_index(network)

    # Per-pathway reachability cache: gene → set of reachable UUIDs
    reachable_cache: dict[str, set[str]] = {}

    # Run all perturbation × key-output predictions
    confusion = defaultdict(int)  # (predicted, expected) → count
    failure_categories: dict[str, int] = defaultdict(int)
    failed_cases: list[dict] = []
    total = 0
    correct = 0
    # A test is "valid" only if both endpoints exist in the network. Drift
    # cases (gene or key output retired between v86 and v96) get tracked
    # separately so we can report accuracy on currently-validatable cases.
    valid_total = 0
    valid_correct = 0

    for gene, direction, col in perturbations:
        gene_uuids = gene_to_uuids.get(gene, [])
        predicted_by_ko: dict[str, int] = {}
        if gene_uuids:
            predicted_by_ko = predict(
                incoming, all_uuids, gene_to_uuids, key_output_uuids, gene, direction,
            )
            if gene not in reachable_cache:
                reachable_cache[gene] = reachable_from(adj, set(gene_uuids))
        else:
            reachable_cache[gene] = set()

        for _, row in curator.iterrows():
            ko = str(row["key_output"])
            ko_uuids = key_output_uuids.get(ko, [])
            if gene_uuids:
                predicted = predicted_by_ko.get(ko, NORMAL) if ko_uuids else NORMAL
            else:
                predicted = NORMAL  # can't perturb what isn't there
            try:
                expected = int(row[col])
            except (ValueError, TypeError):
                continue  # skip malformed cells
            confusion[(predicted, expected)] += 1
            total += 1
            is_valid = bool(gene_uuids) and bool(ko_uuids)
            if is_valid:
                valid_total += 1
            if predicted == expected:
                correct += 1
                if is_valid:
                    valid_correct += 1
            else:
                cat = categorize_failure(
                    predicted, expected, gene_uuids, ko_uuids, reachable_cache[gene],
                )
                failure_categories[cat] += 1
                failed_cases.append({
                    "pathway": pathway_name,
                    "gene": gene,
                    "direction": direction,
                    "key_output": ko,
                    "predicted": predicted,
                    "expected": expected,
                    "category": cat,
                })

    return {
        "status": "ok",
        "name": pathway_name,
        "dbid": pathway_dbid,
        "total": total,
        "correct": correct,
        "accuracy": correct / total if total else 0.0,
        "valid_total": valid_total,
        "valid_correct": valid_correct,
        "valid_accuracy": valid_correct / valid_total if valid_total else 0.0,
        "confusion": dict(confusion),
        "failure_categories": dict(failure_categories),
        "failed_cases": failed_cases,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="output", help="Where regenerated pathway dirs live")
    ap.add_argument("--report", default="/tmp/mpbio_validation.tsv", help="Per-pathway report tsv")
    args = ap.parse_args()

    output_root = Path(args.output_dir)
    pathways = pd.read_csv("/tmp/mpbio_pathways.tsv", sep="\t")
    graph = get_graph()

    rows = []
    overall_confusion: dict[tuple, int] = defaultdict(int)
    overall_failures: dict[str, int] = defaultdict(int)
    all_failed_cases: list[dict] = []
    overall_total = 0
    overall_correct = 0
    overall_valid_total = 0
    overall_valid_correct = 0

    for _, p in pathways.iterrows():
        pid = p["id"]
        name = p["pathway_name"]
        pathway_dir = find_pathway_output_dir(output_root, pid)
        if pathway_dir is None:
            logger.warning(f"No output dir for {name} ({pid}); skipping")
            rows.append({"pathway": name, "status": "no_output", "total": 0, "correct": 0, "accuracy": 0.0})
            continue

        try:
            result = validate_one_pathway(pathway_dir, name, pid, graph)
        except Exception as e:
            logger.error(f"Failed validation for {name}: {e}", exc_info=True)
            rows.append({"pathway": name, "status": f"error: {e}", "total": 0, "correct": 0, "accuracy": 0.0})
            continue

        if result["status"] != "ok":
            rows.append({"pathway": name, "status": result["status"], "total": 0, "correct": 0, "accuracy": 0.0})
            continue

        rows.append({
            "pathway": name,
            "status": "ok",
            "total": result["total"],
            "correct": result["correct"],
            "accuracy": result["accuracy"],
            "valid_total": result["valid_total"],
            "valid_correct": result["valid_correct"],
            "valid_accuracy": result["valid_accuracy"],
        })
        for k, v in result["confusion"].items():
            overall_confusion[k] += v
        for cat, n in result["failure_categories"].items():
            overall_failures[cat] += n
        all_failed_cases.extend(result["failed_cases"])
        overall_total += result["total"]
        overall_correct += result["correct"]
        overall_valid_total += result["valid_total"]
        overall_valid_correct += result["valid_correct"]

    df = pd.DataFrame(rows)
    df.to_csv(args.report, sep="\t", index=False)
    print(df.to_string(index=False))
    print()
    print(f"=== OVERALL: {overall_correct}/{overall_total} = "
          f"{(overall_correct / overall_total * 100 if overall_total else 0):.2f}% raw accuracy ===")
    drift_skipped = overall_total - overall_valid_total
    print(f"=== VALID-ONLY (drift removed): {overall_valid_correct}/{overall_valid_total} = "
          f"{(overall_valid_correct / overall_valid_total * 100 if overall_valid_total else 0):.2f}% ===")
    print(f"    ({drift_skipped} cases skipped because gene or key-output absent in current network)")
    print(f"Confusion (predicted, expected) → count:")
    for (pred, exp), n in sorted(overall_confusion.items()):
        print(f"  pred={pred}, exp={exp}: {n}")
    print()
    print(f"Failure categorization (network bug vs propagator limit):")
    fail_total = sum(overall_failures.values())
    for cat in sorted(overall_failures, key=overall_failures.get, reverse=True):
        n = overall_failures[cat]
        print(f"  {cat}: {n} ({n / fail_total * 100:.1f}%)")
    print(f"Per-pathway report saved to {args.report}")
    failures_path = args.report.replace(".tsv", "_failures.tsv")
    pd.DataFrame(all_failed_cases).to_csv(failures_path, sep="\t", index=False)
    print(f"Per-case failure detail saved to {failures_path}")


if __name__ == "__main__":
    main()
