"""LNG_SUMMATION_REGULATIONS (deltasignal specs/102): regulations a reaction's summation states.

Three Release 97 reactions in the cell cycle make a product no reaction uses, while
their summations say what it does (specs/101 T1, T3). Each rule below writes that
product as a regulator of the reactions the summation names, chosen by structure,
and only inside a pathway that holds both the source reaction and the target.

- E1: RB1:E2F (R-HSA-9018017's product) inhibits every reaction positively
  regulated by an E2F1/2/3 entity. "RB1 binding inhibits transcriptional
  activity of E2F1/2/3:TFDP1/2 complexes."
- E2: AKT's p-T-CDKN1A/B (R-HSA-198613) inhibits every binding of CDKN1A/B to a
  CDK2/4/6 complex. "retention in the cytoplasm, segregating these CDK
  inhibitors from cyclin-CDK complexes."
- E3: CDK4/6:CCND:p-Y-Cip/Kip (R-HSA-8942607) activates every reaction catalysed
  by CDK4/6 with no other active unit. "increased catalytic activity of CDK4 or
  CDK6."

This is a curation supplement to the release, not the release. On by default since
specs/102 (Adam, 2026-10-24); `LNG_SUMMATION_REGULATIONS=0` builds the release alone.
"""
from collections import Counter
from typing import Dict, List, Sequence

STATS: Counter = Counter()

_HOLDS = "-[:hasComponent|hasMember|hasCandidate*0..]->(:EntityWithAccessionedSequence)-[:referenceEntity]->"

# (id, source reaction, regulator, sign, Cypher returning target stIds among $rids)
RULES = [
    ("E1", "R-HSA-9018017", "R-HSA-68644", "neg",
     "UNWIND $rids AS rid MATCH (r:ReactionLikeEvent {stId: rid})"
     "-[:regulatedBy]->(:PositiveRegulation)-[:regulator]->(reg) "
     f"MATCH (reg){_HOLDS}(g:ReferenceGeneProduct) WHERE g.geneName[0] IN ['E2F1', 'E2F2', 'E2F3'] "
     "RETURN DISTINCT r.stId AS target"),
    ("E2", "R-HSA-198613", "R-HSA-198605", "neg",
     "UNWIND $rids AS rid MATCH (r:ReactionLikeEvent {stId: rid})-[:input]->(a) "
     f"MATCH (a){_HOLDS}(ga:ReferenceGeneProduct) WHERE ga.geneName[0] IN ['CDKN1A', 'CDKN1B'] "
     "MATCH (r)-[:input]->(b) WHERE b <> a "
     f"MATCH (b){_HOLDS}(gb:ReferenceGeneProduct) WHERE gb.geneName[0] IN ['CDK2', 'CDK4', 'CDK6'] "
     "RETURN DISTINCT r.stId AS target"),
    ("E3", "R-HSA-8942607", "R-HSA-8942622", "pos",
     "UNWIND $rids AS rid MATCH (r:ReactionLikeEvent {stId: rid})"
     "-[:catalystActivity]->(ca)-[:physicalEntity]->(c) "
     f"MATCH (c){_HOLDS}(g:ReferenceGeneProduct) WHERE g.geneName[0] IN ['CDK4', 'CDK6'] "
     "OPTIONAL MATCH (ca)-[:activeUnit]->(u) "
     "OPTIONAL MATCH (u)" + _HOLDS + "(gu:ReferenceGeneProduct) WHERE gu.geneName[0] IN ['CDK4', 'CDK6'] "
     "WITH r, u, gu WHERE u IS NULL OR gu IS NOT NULL "
     "RETURN DISTINCT r.stId AS target"),
]


def summation_regulations(graph, reaction_ids: Sequence[str]) -> Dict[str, Dict[str, List[str]]]:
    """target reaction stId -> {"pos": [regulator stIds], "neg": [...]} for the
    pathway whose reactions are `reaction_ids`."""
    STATS.clear()
    rids = sorted(set(str(r) for r in reaction_ids))
    held = set(rids)
    out: Dict[str, Dict[str, List[str]]] = {}
    for rule, source, regulator, sign, cypher in RULES:
        if source not in held:
            continue
        for row in graph.run(cypher, rids=rids).data():
            t = row["target"]
            if t == source:
                continue
            lst = out.setdefault(t, {"pos": [], "neg": []})[sign]
            if regulator not in lst:
                lst.append(regulator)
                STATS[rule] += 1
    for d in out.values():
        d["pos"].sort(); d["neg"].sort()
    return out
