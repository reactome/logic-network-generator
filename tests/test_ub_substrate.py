"""LNG_UB_SUBSTRATE_CHANGED (deltasignal specs/057): a ubiquitination depletes
the partner it changes, not the partners it carries through unchanged.

Fixture: "CBL monoubiquitinates activated MET" (R-HSA-8875183) in miniature.
The input complex holds p-MET, GRB2 and p-CBL; the output returns the same
complex with MET ubiquitinated (a new stId) and GRB2 and p-CBL unchanged.
"""
import pandas as pd
import pytest

import src.logic_network_generator as lng
import src.neo4j_connector as nc

UB = "R-HSA-113595"
RX = "R-HSA-RX"
IN_LEAVES = [("R-HSA-pMET", ["MET"]), ("R-HSA-GRB2", ["GRB2"]), ("R-HSA-pCBL", ["CBL"])]
OUT_LEAVES = [("R-HSA-ubMET", ["MET"]), ("R-HSA-GRB2", ["GRB2"]), ("R-HSA-pCBL", ["CBL"])]


class _Result:
    def __init__(self, rows):
        self._rows = rows

    def data(self):
        return self._rows


class _Graph:
    def run(self, query, **kw):
        if "RETURN DISTINCT id AS rxn" in query:
            return _Result([{"rxn": RX}])
        if "-[:input]->(inp:PhysicalEntity)" in query:
            return _Result([{"rxn": RX, "leaf_stid": s, "genes": g} for s, g in IN_LEAVES])
        if "-[:output]->(o:PhysicalEntity)" in query:
            return _Result([{"rxn": RX, "leaf_stid": s, "genes": g} for s, g in OUT_LEAVES])
        raise AssertionError(query)


def _e(a, b, et, pn="pos", ao="and"):
    return {"source_id": a, "target_id": b, "pos_neg": pn, "and_or": ao, "edge_type": et, "stoichiometry": 1}


def _depletions(monkeypatch):
    monkeypatch.setattr(nc, "get_graph", lambda: _Graph())
    monkeypatch.setattr(lng, "_cofactor_stids", lambda: frozenset())
    # The catalyst complex catalyses rx; free p-MET, GRB2 and p-CBL are nodes
    # elsewhere in the network.
    data = [_e("u-cat", "rx", "catalyst"), _e("u-cplx", "rx", "input"), _e("u-ub", "rx", "input"),
            _e("rx", "u-out", "output"),
            _e("u-met", "rx2", "input"), _e("u-grb2", "rx2", "input"), _e("u-cbl", "rx2", "input")]
    r2u = {"u-cat": "R-HSA-CAT", "u-cplx": "R-HSA-CPLX", "u-ub": UB, "u-out": "R-HSA-OUT",
           "u-met": "R-HSA-pMET", "u-grb2": "R-HSA-GRB2", "u-cbl": "R-HSA-pCBL"}
    cmap = pd.DataFrame([{"reaction_uuid": "rx", "reaction_id": RX, "entity_id": "R-HSA-CAT"}])
    lng._emit_substrate_depletion_edges(data, r2u, cmap)
    return {(e["source_id"], e["target_id"]) for e in data if e["edge_type"] == "depletion"}


def test_off_is_the_old_rule_every_partner_is_a_substrate(monkeypatch):
    monkeypatch.delenv("LNG_UB_SUBSTRATE_CHANGED", raising=False)
    assert _depletions(monkeypatch) == {("u-cat", "u-met"), ("u-cat", "u-grb2"), ("u-cat", "u-cbl")}


def test_on_only_the_changed_partner_is_depleted(monkeypatch):
    monkeypatch.setenv("LNG_UB_SUBSTRATE_CHANGED", "1")
    assert _depletions(monkeypatch) == {("u-cat", "u-met")}


def test_the_flag_is_off_by_default():
    from src.env_flags import BOOL_FLAGS
    assert BOOL_FLAGS["LNG_UB_SUBSTRATE_CHANGED"] is False


@pytest.mark.parametrize("changed_only, expect", [
    (False, {"R-HSA-pMET", "R-HSA-GRB2", "R-HSA-pCBL"}),
    (True, {"R-HSA-pMET"}),
])
def test_substrates(changed_only, expect):
    out_genes = {g for _, gs in OUT_LEAVES for g in gs}
    out_stids = {s for s, _ in OUT_LEAVES}
    genes, stids = lng._ubiquitin_substrates(IN_LEAVES, out_genes, out_stids, changed_only)
    assert stids == expect
    assert genes == {g for s, gs in IN_LEAVES if s in expect for g in gs}


def test_nothing_changed_means_no_substrate():
    # Every input leaf reappears unchanged: no depletion at all under the rule
    # (and None, not an empty pair, so the reaction is skipped as before).
    out_genes = {"MET", "GRB2"}
    leaves = [("R-HSA-pMET", ["MET"]), ("R-HSA-GRB2", ["GRB2"])]
    assert lng._ubiquitin_substrates(leaves, out_genes, {"R-HSA-pMET", "R-HSA-GRB2"}, True) is None
    assert lng._ubiquitin_substrates(leaves, set(), set(), False) is None


def test_a_leaf_whose_gene_is_not_in_the_output_is_never_a_substrate():
    leaves = [("R-HSA-X", ["X"])]
    assert lng._ubiquitin_substrates(leaves, {"Y"}, {"R-HSA-Y"}, False) is None
    assert lng._ubiquitin_substrates(leaves, {"Y"}, {"R-HSA-Y"}, True) is None
