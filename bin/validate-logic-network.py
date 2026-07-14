#!/usr/bin/env python
"""Validate a generated logic-network output dir against the LinkML schema.

The CSVs are the flat serialization; multivalued slots (member_leaves,
source_sets, chosen_members) are pipe-delimited cells. This loads them back
into the native LinkML data model (real arrays) and validates the assembled
``LogicNetwork`` instance against schema/logic_network.linkml.yaml.

Usage:
    poetry run python bin/validate-logic-network.py <pathway_output_dir>
"""
import sys
from pathlib import Path

import pandas as pd
from linkml.validator import validate

SCHEMA = Path(__file__).resolve().parent.parent / "schema" / "logic_network.linkml.yaml"
MULTIVALUED = ("member_leaves", "source_sets", "chosen_members")


def _split(v):
    if v is None or (isinstance(v, float) and pd.isna(v)) or v == "":
        return []
    return [p for p in str(v).split("|") if p]


def _clean(v):
    return None if (v is None or (isinstance(v, float) and pd.isna(v)) or v == "") else v


def load_instance(d: Path) -> dict:
    nodes = []
    for r in pd.read_csv(d / "nodes.csv", dtype=str, keep_default_na=False).to_dict("records"):
        node = {k: _clean(v) for k, v in r.items() if k not in MULTIVALUED}
        for k in MULTIVALUED:
            node[k] = _split(r.get(k))
        nodes.append({k: v for k, v in node.items() if v is not None or k in MULTIVALUED})

    contexts = pd.read_csv(d / "node_reaction_context.csv", dtype=str,
                           keep_default_na=False).to_dict("records")

    edges = []
    for r in pd.read_csv(d / "logic_network.csv", keep_default_na=False).to_dict("records"):
        e = {k: _clean(v) for k, v in r.items()}
        if e.get("stoichiometry") not in (None, ""):
            e["stoichiometry"] = int(float(e["stoichiometry"]))
        edges.append({k: v for k, v in e.items() if v is not None})

    return {"pathway_id": d.name, "nodes": nodes,
            "node_reaction_contexts": contexts, "edges": edges}


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    d = Path(sys.argv[1])
    report = validate(load_instance(d), str(SCHEMA), "LogicNetwork")
    if not report.results:
        print(f"VALID: {d} conforms to {SCHEMA.name}")
        return 0
    print(f"INVALID: {len(report.results)} problem(s) in {d}")
    for res in report.results[:25]:
        print(f"  [{res.severity}] {res.message}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
