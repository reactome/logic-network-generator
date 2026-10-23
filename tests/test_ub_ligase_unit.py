"""deltasignal specs/089: a ubiquitination does not deplete its own ligase
(LNG_UB_LIGASE_UNIT)."""
from typing import Any, Dict

import pandas as pd
import pytest

import src.logic_network_generator as lng
import src.neo4j_connector as nc

UB = "R-HSA-113595"
RX = "R-HSA-1363331"          # SCF(Skp2) ubiquitinates p130 and is released
CKS1B, RBL2, UB_RBL2 = "R-HSA-CKS1B", "R-HSA-RBL2", "R-HSA-UB-RBL2"


def _e(a, b, et, pn="pos", ao="and"):
    return {"source_id": a, "target_id": b, "pos_neg": pn, "and_or": ao, "edge_type": et,
            "stoichiometry": 1}


# --- the substrate rule on its own ------------------------------------------

def test_a_ligase_leaf_carried_unchanged_is_not_the_substrate():
    genes, stids = lng._ub_substrates([(CKS1B, ["CKS1B"]), (RBL2, ["RBL2"])],
                                      {"CKS1B", "RBL2"}, {CKS1B, UB_RBL2}, {CKS1B})
    assert stids == {RBL2} and genes == {"RBL2"}


def test_without_ligase_leaves_every_shared_gene_is_the_substrate():
    _, stids = lng._ub_substrates([(CKS1B, ["CKS1B"]), (RBL2, ["RBL2"])],
                                  {"CKS1B", "RBL2"}, {CKS1B, UB_RBL2}, set())
    assert stids == {CKS1B, RBL2}


def test_an_autoubiquitinated_ligase_stays_the_substrate():
    # BRAP autoubiquitinates: BRAP is the active unit, but its output form is
    # Ub-BRAP, another stId, so it is modified and still depleted.
    _, stids = lng._ub_substrates([("R-HSA-BRAP", ["BRAP"])], {"BRAP"},
                                  {"R-HSA-UB-BRAP"}, {"R-HSA-BRAP"})
    assert stids == {"R-HSA-BRAP"}


def test_the_substrates_partners_are_kept():
    # CBL ubiquitinates EGFR: CBL (active unit) is dropped, the receptor complex
    # partners are still depleted (specs/078 lost the CBL class by dropping them).
    ins = [("R-HSA-CBL", ["CBL"]), ("R-HSA-EGFR", ["EGFR"]), ("R-HSA-GRB2", ["GRB2"])]
    _, stids = lng._ub_substrates(ins, {"CBL", "EGFR", "GRB2"},
                                  {"R-HSA-CBL", "R-HSA-EGFR", "R-HSA-GRB2"}, {"R-HSA-CBL"})
    assert stids == {"R-HSA-EGFR", "R-HSA-GRB2"}


# --- the whole pass, against a stub graph ------------------------------------

class _Result:
    def __init__(self, rows): self.rows = rows
    def data(self): return self.rows


class _Graph:
    """Answers the ubiquitin pass's queries for one reaction: catalyst C
    (p-RBL2:SCF(Skp2):Cks1) + Ub -> Ub-p-RBL2 + the released ligase, no activeUnit."""
    def run(self, q: str, **kw: Any) -> _Result:
        if "activeUnit" in q:
            return _Result([])
        if "<-[:hasComponent]-(:Complex)<-[:input]-" in q:
            return _Result([{"rxn": RX, "leaf_stid": CKS1B}])
        if "-[:output]->" in q:
            return _Result([{"rxn": RX, "leaf_stid": CKS1B, "genes": ["CKS1B"]},
                            {"rxn": RX, "leaf_stid": UB_RBL2, "genes": ["RBL2"]}])
        if "-[:input]->(inp" in q:
            return _Result([{"rxn": RX, "leaf_stid": CKS1B, "genes": ["CKS1B"]},
                            {"rxn": RX, "leaf_stid": RBL2, "genes": ["RBL2"]}])
        if "ub.stId IN $ub_stids" in q:
            return _Result([{"rxn": RX}])
        raise AssertionError(q)


@pytest.fixture
def pass_(monkeypatch):
    monkeypatch.setattr(nc, "get_graph", lambda: _Graph())
    monkeypatch.setattr(lng, "_cofactor_stids", lambda: frozenset())
    monkeypatch.setenv("LNG_VARIANT_NODES", "0")

    def run() -> set:
        data = [_e("u-c", "rx", "catalyst"), _e("u-c", "rx", "input"), _e("rx", "u-out", "output"),
                _e("u-cks", "x1", "input"), _e("u-rbl2", "x2", "input")]
        r2u: Dict[str, str] = {"u-c": "R-HSA-1363326", "u-out": UB_RBL2,
                               "u-cks": CKS1B, "u-rbl2": RBL2}
        cmap = pd.DataFrame([{"reaction_uuid": "rx", "reaction_id": RX}])
        lng._emit_substrate_depletion_edges(data, r2u, cmap)
        return {(e["source_id"], e["target_id"]) for e in data if e["edge_type"] == "depletion"}
    return run


def test_off_the_ligase_drains_its_own_subunit(pass_, monkeypatch):
    monkeypatch.setenv("LNG_UB_LIGASE_UNIT", "0")
    assert pass_() == {("u-c", "u-cks"), ("u-c", "u-rbl2")}


def test_on_only_the_substrate_is_drained(pass_, monkeypatch):
    monkeypatch.setenv("LNG_UB_LIGASE_UNIT", "1")
    assert pass_() == {("u-c", "u-rbl2")}


def test_the_rule_is_the_default(pass_, monkeypatch):
    monkeypatch.delenv("LNG_UB_LIGASE_UNIT", raising=False)
    assert lng.env_flag("LNG_UB_LIGASE_UNIT") is True
    assert pass_() == {("u-c", "u-rbl2")}


def test_rule_b_applies_only_where_no_active_unit_is_curated(monkeypatch):
    # review 2026-10-22: the stub above answers activeUnit with nothing, so rule A
    # and its precedence over rule B were untested. R1 has an active unit (the
    # ligase L); its release query would also name a partner P, which must stay.
    class G:
        def run(self, q, **kw):
            if "activeUnit" in q:
                return _Result([{"rxn": "R1", "leaf_stid": "L"}])
            if "<-[:hasComponent]-(:Complex)<-[:input]-" in q:
                assert kw["ids"] == ["R2"]          # R1 is not asked again
                rows = [{"rxn": "R1", "leaf_stid": "P"}, {"rxn": "R2", "leaf_stid": "Q"}]
                return _Result([r for r in rows if r["rxn"] in kw["ids"]])
            raise AssertionError(q)
    monkeypatch.setattr(nc, "get_graph", lambda: G())
    assert lng._ub_ligase_leaves(["R1", "R2"]) == {"R1": {"L"}, "R2": {"Q"}}
