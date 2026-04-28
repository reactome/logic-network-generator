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


def propagate(
    network: pd.DataFrame,
    pinned: dict[str, int],
) -> dict[str, int]:
    """Run Boolean propagation until fixed point or MAX_ITERATIONS."""
    # Build incoming-edge index per target uuid
    incoming = defaultdict(list)
    all_uuids: set[str] = set()
    for _, row in network.iterrows():
        target = row["target_id"]
        source = row["source_id"]
        all_uuids.add(target)
        all_uuids.add(source)
        incoming[target].append({
            "source": source,
            "edge_type": row["edge_type"],
            "pos_neg": row["pos_neg"],
        })

    state = {u: pinned.get(u, NORMAL) for u in all_uuids}

    for _ in range(MAX_ITERATIONS):
        new_state = dict(state)
        for uuid_, edges in incoming.items():
            if uuid_ in pinned:
                continue  # perturbed nodes stay pinned
            and_contribs: list[int] = []
            or_contribs: list[int] = []
            for e in edges:
                src = state.get(e["source"], NORMAL)
                etype = e["edge_type"]
                if etype == "regulator" and e["pos_neg"] == "neg":
                    and_contribs.append(invert(src))
                elif etype in {"input", "catalyst", "regulator", "assembly"}:
                    and_contribs.append(src)
                elif etype in {"output", "dissociation"}:
                    or_contribs.append(src)

            if and_contribs and or_contribs:
                new_state[uuid_] = max(min(and_contribs), max(or_contribs))
            elif and_contribs:
                new_state[uuid_] = min(and_contribs)
            elif or_contribs:
                new_state[uuid_] = max(or_contribs)
            # else: no edges, keep prior state
        if new_state == state:
            break
        state = new_state
    return state


def predict(network, gene_to_uuids, key_output_uuids, perturbation_gene, direction):
    """Run one perturbation; return predicted states for each key-output uuid."""
    pinned: dict[str, int] = {}
    for u in gene_to_uuids.get(perturbation_gene, []):
        pinned[u] = direction
    final = propagate(network, pinned)
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


def validate_one_pathway(pathway_dir: Path, pathway_name: str, pathway_dbid: str, graph) -> dict:
    """Validate one pathway, return per-pathway metrics."""
    curator_path = MP_BIOPATH_DIR / "reactome_curator_predictions" / f"{pathway_name}_reactome_curator_results.tsv"
    if not curator_path.exists():
        return {"status": "no_curator_file", "name": pathway_name}

    curator = pd.read_csv(curator_path, sep="\t")
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

    # Run all perturbation × key-output predictions
    confusion = defaultdict(int)  # (predicted, expected) → count
    total = 0
    correct = 0
    skipped_missing_uuid = 0
    skipped_missing_gene = 0

    for gene, direction, col in perturbations:
        if gene not in gene_to_uuids or not gene_to_uuids[gene]:
            skipped_missing_gene += len(curator)
            continue
        predicted_by_ko = predict(
            network, gene_to_uuids, key_output_uuids, gene, direction
        )
        for _, row in curator.iterrows():
            ko = str(row["key_output"])
            if not key_output_uuids.get(ko):
                skipped_missing_uuid += 1
                continue
            predicted = predicted_by_ko[ko]
            expected = int(row[col])
            confusion[(predicted, expected)] += 1
            total += 1
            if predicted == expected:
                correct += 1

    return {
        "status": "ok",
        "name": pathway_name,
        "dbid": pathway_dbid,
        "total": total,
        "correct": correct,
        "accuracy": correct / total if total else 0.0,
        "confusion": dict(confusion),
        "skipped_missing_gene": skipped_missing_gene,
        "skipped_missing_uuid": skipped_missing_uuid,
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
    overall_total = 0
    overall_correct = 0

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
        })
        for k, v in result["confusion"].items():
            overall_confusion[k] += v
        overall_total += result["total"]
        overall_correct += result["correct"]

    df = pd.DataFrame(rows)
    df.to_csv(args.report, sep="\t", index=False)
    print(df.to_string(index=False))
    print()
    print(f"=== OVERALL: {overall_correct}/{overall_total} = "
          f"{(overall_correct / overall_total * 100 if overall_total else 0):.2f}% accuracy ===")
    print(f"Confusion (predicted, expected) → count:")
    for (pred, exp), n in sorted(overall_confusion.items()):
        print(f"  pred={pred}, exp={exp}: {n}")
    print(f"Report saved to {args.report}")


if __name__ == "__main__":
    main()
