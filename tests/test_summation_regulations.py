"""deltasignal specs/102: regulations a reaction's summation states
(LNG_SUMMATION_REGULATIONS)."""
import pytest

import src.summation_regulations as sr


class _FakeGraph:
    def __init__(self, targets):
        self.targets = targets

    def run(self, cypher, rids):
        rule = next(r for r in sr.RULES if r[4] == cypher)[0]
        rows = [{"target": t} for t in self.targets.get(rule, []) if t in rids]

        class _R:
            def data(self):
                return rows
        return _R()


def test_a_rule_applies_only_where_its_source_reaction_is():
    g = _FakeGraph({"E1": ["R-HSA-TK1"], "E2": ["R-HSA-BIND"]})
    # Mitotic G1 holds RB1 binding E2F (E1's source) but not AKT's p27 phosphorylation (E2's)
    regs = sr.summation_regulations(g, ["R-HSA-9018017", "R-HSA-TK1", "R-HSA-BIND"])
    assert regs == {"R-HSA-TK1": {"pos": [], "neg": ["R-HSA-68644"]}}
    assert sr.STATS == {"E1": 1}


def test_signs_and_the_source_is_never_its_own_target():
    g = _FakeGraph({"E2": ["R-HSA-198613", "R-HSA-BIND"], "E3": ["R-HSA-69227"]})
    regs = sr.summation_regulations(g, ["R-HSA-198613", "R-HSA-8942607", "R-HSA-BIND", "R-HSA-69227"])
    assert regs == {"R-HSA-BIND": {"pos": [], "neg": ["R-HSA-198605"]},
                    "R-HSA-69227": {"pos": ["R-HSA-8942622"], "neg": []}}


def test_flag_off_by_default(monkeypatch):
    from src.env_flags import env_flag
    monkeypatch.delenv("LNG_SUMMATION_REGULATIONS", raising=False)
    assert env_flag("LNG_SUMMATION_REGULATIONS") is False


@pytest.mark.database
def test_release_targets_in_mitotic_g1():
    """The registered targets (specs/102 spec.md), read from the database."""
    from src.neo4j_connector import get_graph
    graph = get_graph()
    rids = [r["s"] for r in graph.run(
        "MATCH (:Pathway {stId: 'R-HSA-453279'})-[:hasEvent*]->(r:ReactionLikeEvent) "
        "RETURN DISTINCT r.stId AS s").data()]
    regs = sr.summation_regulations(graph, rids)
    e1 = [t for t, d in regs.items() if "R-HSA-68644" in d["neg"]]
    assert len(e1) == 15
    assert sorted(t for t, d in regs.items() if "R-HSA-198605" in d["neg"]) == ["R-HSA-69562", "R-HSA-8941915"]
    assert sorted(t for t, d in regs.items() if "R-HSA-8942622" in d["pos"]) == [
        "R-HSA-1226094", "R-HSA-1226095", "R-HSA-69227"]
