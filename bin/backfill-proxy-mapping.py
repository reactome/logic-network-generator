#!/usr/bin/env python3
"""Backfill entity_reaction_proxy_mapping.csv for already-generated pathways.

Regenerating a network from scratch is expensive; this reuses the existing
logic_network.csv + stid_to_uuid_mapping.csv and only adds the new proxy file
(which just needs Neo4j + those two artifacts). Going forward the generator
emits the file itself; this is a one-time catch-up for the existing catalog.
"""
import sys
from pathlib import Path

import pandas as pd

from src.neo4j_connector import prefetch_entity_data
from src.logic_network_generator import export_entity_reaction_proxy_mapping

OUTPUT = Path("output")


def backfill(pathway_dir: Path) -> None:
    name = pathway_dir.name
    if not name.endswith(tuple(f"R-HSA-{n}" for n in [""])) and "R-HSA-" not in name:
        return
    pid = "R-HSA-" + name.rsplit("R-HSA-", 1)[-1]
    net_f = pathway_dir / "logic_network.csv"
    map_f = pathway_dir / "stid_to_uuid_mapping.csv"
    rc_f = pathway_dir / "cache" / "reaction_connections.csv"
    if not (net_f.exists() and map_f.exists()):
        print(f"  skip {name}: missing artifacts")
        return

    network = pd.read_csv(net_f)
    mapping = pd.read_csv(map_f)
    reactome_id_to_uuid = dict(zip(mapping["uuid"].astype(str),
                                   mapping["stable_id"].astype(str)))
    # Approximate reaction_id_map from the mapping (reaction stIds carry UUIDs
    # too); only reaction stIds get queried, so entity rows are harmless.
    rid_map = mapping.rename(columns={"uuid": "uid", "stable_id": "reactome_id"})

    # Warm caches with the pathway's reactions for fast terminal-component walks.
    if rc_f.exists():
        rc = pd.read_csv(rc_f, dtype=str)
        rxn_ids = set(rc["preceding_reaction_id"].dropna()) | \
                  set(rc["following_reaction_id"].dropna())
        if rxn_ids:
            prefetch_entity_data(list(rxn_ids))

    out_f = pathway_dir / "entity_reaction_proxy_mapping.csv"
    export_entity_reaction_proxy_mapping(
        network, rid_map, reactome_id_to_uuid, pid, str(out_f))
    df = pd.read_csv(out_f)
    print(f"  {name}: {len(df)} rows, "
          f"{df['entity_stable_id'].nunique() if len(df) else 0} entities")


if __name__ == "__main__":
    targets = sys.argv[1:]
    if targets:
        dirs = [OUTPUT / t for t in targets]
    else:
        dirs = [d for d in sorted(OUTPUT.iterdir()) if d.is_dir()]
    for d in dirs:
        try:
            backfill(d)
        except Exception as e:
            print(f"  ERROR {d.name}: {e}")
