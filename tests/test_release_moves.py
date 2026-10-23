"""deltasignal specs/093: carrier-loop releases that move the free form across
compartments (a receptor recycled to the membrane) are marked in carrier_loops.csv."""
import src.logic_network_generator as lng
import src.neo4j_connector as nc


class _Result:
    def __init__(self, rows): self.rows = rows
    def data(self): return self.rows


class _Graph:
    def run(self, q, **kw):
        assert "compartment" in q
        rows = {("R-MET", "R-RAB4"): (["PM"], ["RE"]),          # recycled: moves
                ("R-MET", "R-STAT3"): (["PM"], ["PM", "CYT"]),  # released in place
                ("R-X", "R-Y"): ([], ["CYT"])}                  # no compartment: not marked
        return _Result([{"f": f, "r": r, "cf": rows[(f, r)][0], "ci": rows[(f, r)][1]}
                        for f, r in kw["p"] if (f, r) in rows])


def test_only_a_release_across_compartments_moves(monkeypatch):
    monkeypatch.setattr(nc, "get_graph", lambda: _Graph())
    assert lng._release_moves({("R-MET", "R-RAB4"), ("R-MET", "R-STAT3"), ("R-X", "R-Y")}) == {("R-MET", "R-RAB4")}


def test_no_pairs_asks_nothing(monkeypatch):
    monkeypatch.setattr(nc, "get_graph", lambda: (_ for _ in ()).throw(AssertionError("queried")))
    assert lng._release_moves(set()) == set()
