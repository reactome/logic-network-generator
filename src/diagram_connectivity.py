"""Reaction connectivity from Reactome diagram JSON.

The generator connects two reactions only when the curator annotated a
``precedingEvent`` between them. Older pathways under-annotate ``precedingEvent``,
so reactions where A's product is literally B's substrate are left disconnected.

The pathway **diagram** is the curator's drawn connectivity: two reactions are
linked when they share an entity **glyph** (A's output glyph == B's input/catalyst
glyph). We extract those (producer, consumer) reaction pairs and feed them into
``reaction_connections`` alongside ``precedingEvent`` — Phase 2 then merges the
shared product into one node, giving ``A -> product -> B``.

Why the diagram beats raw Neo4j input/output matching: curators draw cofactors
(ATP/ADP/H2O/etc.) as **separate per-reaction glyphs**, so shared-glyph
connectivity excludes cofactors for free — no hub/threshold tuning.

See reactome/logic-network-generator#39.
"""
import json
import os
from collections import defaultdict
from pathlib import Path
from typing import Set, Tuple

import pandas as pd

from src.argument_parser import logger


def _diagram_dir() -> Path:
    return Path(os.environ.get("LNG_DIAGRAM_DIR", os.path.expanduser("~/reactome-diagrams/97")))


def _has_diagram(stid: str) -> bool:
    d = _diagram_dir()
    return (d / f"{stid}.json").exists() and (d / f"{stid}.graph.json").exists()


def _covering_diagram_stid(pathway_id: str) -> str:
    """stId of the diagram that renders this pathway.

    We generate mostly top-level pathways (which have their own diagram), but the
    TEST pathways are often sub-pathways that inherit an ancestor's diagram. Use
    the pathway's own diagram if present, else the nearest diagrammed ancestor
    (walking up ``hasEvent`` in Neo4j). Returns "" if none found.
    """
    if _has_diagram(pathway_id):
        return pathway_id
    from src.neo4j_connector import get_graph
    ancestors = get_graph().run(
        """MATCH path=(anc:Pathway)-[:hasEvent*]->(p:Pathway {stId:$pid})
           RETURN anc.stId AS stid, length(path) AS d ORDER BY d""",
        pid=pathway_id,
    ).data()
    for a in ancestors:
        if _has_diagram(a["stid"]):
            return a["stid"]
    return ""


def _pathway_reaction_stids(pathway_id: str) -> Set[str]:
    """All ReactionLikeEvent stIds contained in the pathway (to isolate it within
    an ancestor diagram)."""
    from src.neo4j_connector import get_graph
    rows = get_graph().run(
        "MATCH (p:Pathway {stId:$pid})-[:hasEvent*]->(r:ReactionLikeEvent) "
        "RETURN collect(DISTINCT r.stId) AS r",
        pid=pathway_id,
    ).evaluate()
    return set(rows or [])


def diagram_shared_product_pairs(pathway_id: str) -> Set[Tuple[str, str]]:
    """(producer_stId, consumer_stId) reaction pairs drawn as connected in the diagram.

    A pair is emitted when a producer reaction's OUTPUT glyph is the same glyph a
    consumer reaction takes as INPUT or CATALYST. Reaction identity is the
    Reactome stId (mapped from the diagram's numeric reactomeId via the
    companion ``.graph.json``). The diagram used is the pathway's own, else the
    nearest diagrammed ancestor; pairs are then **restricted to reactions that
    belong to this pathway** so an ancestor diagram only contributes the target
    pathway's internal connectivity. Returns an empty set if no diagram covers it.
    """
    ddir = _diagram_dir()
    diagram_stid = _covering_diagram_stid(pathway_id)
    if not diagram_stid:
        logger.info(f"No diagram (own or ancestor) covers {pathway_id}; skipping diagram connectivity")
        return set()
    if diagram_stid != pathway_id:
        logger.info(f"{pathway_id} has no own diagram; using ancestor diagram {diagram_stid}")

    layout = json.loads((ddir / f"{diagram_stid}.json").read_text())
    graph = json.loads((ddir / f"{diagram_stid}.graph.json").read_text())
    # Restrict to this pathway's reactions (isolate it within the ancestor diagram).
    own_reactions = _pathway_reaction_stids(pathway_id)

    # reaction dbId -> stId (graph.json edges carry both)
    dbid_to_stid = {e["dbId"]: e["stId"] for e in graph.get("edges", []) if e.get("stId")}

    # For each reaction glyph (diagram 'edge'): its output glyph ids, and its
    # input+catalyst glyph ids. reactomeId on a diagram edge is the reaction dbId.
    glyph_out_rxns = defaultdict(set)   # glyph_id -> {reaction_dbId producing it}
    glyph_in_rxns = defaultdict(set)    # glyph_id -> {reaction_dbId consuming it}
    for e in layout.get("edges", []):
        rdb = e.get("reactomeId")
        if rdb is None:
            continue
        for x in e.get("outputs", []):
            glyph_out_rxns[x["id"]].add(rdb)
        for role in ("inputs", "catalysts"):
            for x in e.get(role, []):
                glyph_in_rxns[x["id"]].add(rdb)

    pairs: Set[Tuple[str, str]] = set()
    for glyph_id in set(glyph_out_rxns) & set(glyph_in_rxns):
        for producer in glyph_out_rxns[glyph_id]:
            for consumer in glyph_in_rxns[glyph_id]:
                if producer == consumer:
                    continue
                p_st = dbid_to_stid.get(producer)
                c_st = dbid_to_stid.get(consumer)
                if not p_st or not c_st:
                    continue
                # Keep only pairs whose BOTH reactions belong to this pathway
                # (so an ancestor diagram contributes only this pathway's flow).
                if own_reactions and (p_st not in own_reactions or c_st not in own_reactions):
                    continue
                pairs.add((p_st, c_st))
    return pairs


def augment_reaction_connections(pathway_id: str,
                                 reaction_connections: pd.DataFrame) -> pd.DataFrame:
    """Union diagram-drawn product->substrate pairs into reaction_connections.

    Adds only pairs not already present (as preceding->following). Tagged
    event_status='Diagram Shared Product' for traceability. No-op (returns input
    unchanged) when disabled or no diagram is available.
    """
    if os.environ.get("LNG_DIAGRAM_CONNECTIVITY", "1") == "0":
        return reaction_connections

    pairs = diagram_shared_product_pairs(pathway_id)
    if not pairs:
        return reaction_connections

    existing = set(
        zip(reaction_connections.get("preceding_reaction_id", pd.Series(dtype=str)),
            reaction_connections.get("following_reaction_id", pd.Series(dtype=str)))
    )
    new_rows = [
        {"preceding_reaction_id": p, "following_reaction_id": c,
         "event_status": "Diagram Shared Product"}
        for (p, c) in pairs if (p, c) not in existing
    ]
    if not new_rows:
        logger.info(f"Diagram connectivity for {pathway_id}: all {len(pairs)} pairs already linked")
        return reaction_connections

    logger.info(
        f"Diagram connectivity for {pathway_id}: +{len(new_rows)} product->substrate "
        f"pairs not in precedingEvent (of {len(pairs)} drawn)"
    )
    return pd.concat([reaction_connections, pd.DataFrame(new_rows)], ignore_index=True)
