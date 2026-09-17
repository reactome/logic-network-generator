"""Layered connectivity tests: source fact -> exported fact -> structure -> reachability.

The end-to-end benchmark tells us a pathway scores badly. It does not tell us
which layer broke. These four tiers each assert the same real defect at a
different depth, so a failure localises itself:

  Tier 1  what Reactome actually says            (needs Neo4j)
  Tier 2  what containment.csv exports about it  (needs Neo4j + artifacts)
  Tier 3  a structural invariant of the network  (artifacts only)
  Tier 4  the end-to-end reachability it implies (artifacts only)

The worked case is the ISGF3 nuclear-import chain in Interferon alpha/beta
signalling (R-HSA-909733). Reactome connects cytosolic ISGF3 to the nuclear
transcription machinery through COMPOSITION ONLY — there is no curated binding
reaction — so a reaction-only logic network cannot cross it. The generator
bridges the gap with synthetic assembly/dissociation edges, but builds the two
sides from different uuids, so the bridge does not join up and every
perturbation upstream of ISGF3 fails to reach the readout.

Tiers 2-4 are xfail(strict) — they encode the behaviour we want and will flag
loudly the moment a fix makes them pass. Tier 1 must pass today: it guards the
premise the other three rest on.
"""

import os
from pathlib import Path

import pandas as pd
import pytest

# --- the worked case, by stable id -----------------------------------------
PATHWAY = "R-HSA-909733"          # Interferon alpha/beta signaling
ISGF3_CYTOSOL = "R-HSA-909698"    # ISGF3 [cytosol]
ISGF3_KPNA1 = "R-HSA-9710958"     # ISGF3:KPNA1 [cytosol]        (intermediate)
ISGF3_IMPORTIN = "R-HSA-9710965"  # ISGF3:KPNA1:KPNB1 [cytosol]
TRANSLOCATION = "R-HSA-909721"    # Translocation of ISGF3 complex to nucleus
ISGF3_NUCLEUS = "R-HSA-913527"    # ISGF3 [nucleoplasm]
EXPRESSION = "R-HSA-1015702"      # Expression of IFN-induced genes (the readout)


def _catalog_root() -> Path:
    """Where generated pathway bundles live.

    Defaults to output/, the repo convention. LNG_TEST_CATALOG points the
    tiers at an existing catalog build instead of regenerating one.
    """
    return Path(os.environ.get("LNG_TEST_CATALOG", "output"))


def _bundle() -> Path | None:
    root = _catalog_root()
    if not root.exists():
        return None
    for d in root.iterdir():
        if d.is_dir() and d.name.endswith(PATHWAY) and (d / "logic_network.csv").exists():
            return d
    return None


BUNDLE = _bundle()
needs_bundle = pytest.mark.skipif(
    BUNDLE is None,
    reason=f"No generated {PATHWAY} bundle under {_catalog_root()} "
           f"(set LNG_TEST_CATALOG to a catalog build)",
)


def _stid_to_uuids(bundle: Path) -> dict[str, set[str]]:
    """stable id -> the uuids standing for it in this network.

    Requires the columns by name. A positional fallback would silently swap
    them — the file ships as `uuid,stable_id`, the opposite of the reading
    order — and every tier below would then measure nonsense while passing.
    """
    m = pd.read_csv(bundle / "stid_to_uuid_mapping.csv", dtype=str)
    missing = {"stable_id", "uuid"} - set(m.columns)
    assert not missing, f"stid_to_uuid_mapping.csv is missing {missing}"
    out: dict[str, set[str]] = {}
    for stid, uuid in zip(m["stable_id"], m["uuid"]):
        out.setdefault(str(stid), set()).add(str(uuid))
    return out


# --- Tier 1: what Reactome actually says ------------------------------------

@pytest.mark.database
class TestReactomeSourceFacts:
    """The premise the other tiers rest on.

    If Reactome later curates the missing binding reactions, these fail and
    tell us the workaround is no longer needed — rather than the workaround
    quietly staying in place.
    """

    @pytest.fixture(scope="class")
    def graph(self):
        from py2neo import Graph  # type: ignore
        return Graph(
            os.environ.get("NEO4J_URI", "bolt://localhost:7687"),
            auth=(os.environ.get("NEO4J_USER", "neo4j"),
                  os.environ.get("NEO4J_PASSWORD", "neo4j")),
        )

    def _consumers(self, graph, stid):
        return [r["s"] for r in graph.run(
            "MATCH (rx:ReactionLikeEvent)-[:input]->(e {stId:$s}) RETURN rx.stId AS s", s=stid)]

    def _producers(self, graph, stid):
        return [r["s"] for r in graph.run(
            "MATCH (rx:ReactionLikeEvent)-[:output]->(e {stId:$s}) RETURN rx.stId AS s", s=stid)]

    def test_cytosolic_isgf3_is_consumed_by_no_reaction(self, graph):
        """The gap itself: nothing downstream of ISGF3 [cytosol] by reaction."""
        assert self._consumers(graph, ISGF3_CYTOSOL) == [], (
            "ISGF3 [cytosol] now has a consuming reaction — the composition-only "
            "gap this module works around may be curated; re-check the bridge."
        )

    def test_importin_complex_is_produced_by_no_reaction(self, graph):
        assert self._producers(graph, ISGF3_IMPORTIN) == [], (
            "ISGF3:KPNA1:KPNB1 now has a producing reaction; the assembly "
            "bridge may no longer be needed."
        )

    def test_the_chain_is_joined_by_composition(self, graph):
        """ISGF3 -> ISGF3:KPNA1 -> ISGF3:KPNA1:KPNB1, via hasComponent."""
        def parents(stid):
            return [r["s"] for r in graph.run(
                "MATCH (c:Complex)-[:hasComponent]->(e {stId:$s}) RETURN c.stId AS s", s=stid)]

        assert ISGF3_KPNA1 in parents(ISGF3_CYTOSOL)
        assert ISGF3_IMPORTIN in parents(ISGF3_KPNA1)

    def test_translocation_carries_the_complex_into_the_nucleus(self, graph):
        assert TRANSLOCATION in self._consumers(graph, ISGF3_IMPORTIN)
        assert TRANSLOCATION in self._producers(graph, ISGF3_NUCLEUS)


# --- Tier 2: what containment.csv exports about it --------------------------

@pytest.mark.database
@pytest.mark.integration
@needs_bundle
class TestContainmentCompleteness:

    @pytest.fixture(scope="class")
    def containment(self):
        return pd.read_csv(BUNDLE / "containment.csv", dtype=str)

    def test_containment_records_the_leaves(self, containment):
        """Passes today — the leaf-level fact IS exported."""
        rows = containment[containment["stable_id"] == ISGF3_IMPORTIN]
        contained = set(rows["contains_stable_id"])
        assert "R-HSA-879183" in contained, "IRF9 should be inside the importin complex"

    @pytest.mark.xfail(
        strict=True,
        reason="export_containment recurses to the LEAVES, so an intermediate "
               "complex that is itself a node in the network is never recorded "
               "as contained. A consumer cannot learn ISGF3 is inside "
               "ISGF3:KPNA1:KPNB1, only that they share subunits.",
    )
    def test_containment_records_intermediate_complexes(self, containment):
        """Containment must be complete over entities that ARE nodes here.

        If X and Y are both nodes in this network and Reactome says X is a
        (transitive) component of Y, containment.csv must say so. Leaf-sharing
        is not a substitute: two unrelated complexes can share every subunit.
        """
        rows = containment[containment["stable_id"] == ISGF3_IMPORTIN]
        contained = set(rows["contains_stable_id"])
        assert ISGF3_CYTOSOL in contained


# --- Tier 3: a structural invariant of the generated network ----------------

@pytest.mark.integration
@needs_bundle
class TestBoundaryBridgeCoherence:

    @pytest.fixture(scope="class")
    def network(self):
        return pd.read_csv(BUNDLE / "logic_network.csv", dtype=str)

    @pytest.fixture(scope="class")
    def uuids(self):
        return _stid_to_uuids(BUNDLE)

    @pytest.mark.xfail(
        strict=True,
        reason="Dissociation targets a freshly minted per-occurrence sink while "
               "assembly sources the shared functional node, so the two sides of "
               "the same molecule never share a uuid and the bridge does not "
               "join up. See the 'Fresh per-occurrence readout sink' comment in "
               "logic_network_generator.py.",
    )
    def test_assembly_and_dissociation_sides_share_a_uuid(self, network, uuids):
        """A molecule released by one complex must be the one another consumes.

        Any stable id appearing on both sides of the boundary bridge must share
        at least one uuid between them — otherwise signal entering the
        dissociation side can never leave through the assembly side.
        """
        uuid_to_stid = {u: s for s, us in uuids.items() for u in us}
        diss = network[network["edge_type"] == "dissociation"]["target_id"]
        asm = network[network["edge_type"] == "assembly"]["source_id"]
        diss_by_stid: dict[str, set[str]] = {}
        asm_by_stid: dict[str, set[str]] = {}
        for u in diss:
            diss_by_stid.setdefault(uuid_to_stid.get(u, u), set()).add(u)
        for u in asm:
            asm_by_stid.setdefault(uuid_to_stid.get(u, u), set()).add(u)

        shared_both_sides = set(diss_by_stid) & set(asm_by_stid)
        assert shared_both_sides, "no stable id on both sides — nothing to check"
        broken = {s for s in shared_both_sides if not (diss_by_stid[s] & asm_by_stid[s])}
        assert not broken, (
            f"{len(broken)} of {len(shared_both_sides)} molecules appear on both "
            f"sides of the boundary bridge with no shared uuid: {sorted(broken)[:5]}"
        )


# --- Tier 4: the end-to-end reachability it implies -------------------------

@pytest.mark.integration
@needs_bundle
class TestNuclearImportReachability:

    @pytest.fixture(scope="class")
    def forward(self):
        net = pd.read_csv(BUNDLE / "logic_network.csv", dtype=str)
        adj: dict[str, set[str]] = {}
        for src, tgt in zip(net["source_id"], net["target_id"]):
            adj.setdefault(src, set()).add(tgt)
        return adj

    @pytest.fixture(scope="class")
    def uuids(self):
        return _stid_to_uuids(BUNDLE)

    def test_every_entity_in_the_chain_is_a_node_here(self, uuids):
        """Guards the xfail below.

        `uuids[X]` on an absent stable id raises KeyError, and an error counts
        as an expected failure just like an assertion does — so the xfail would
        look satisfied for entirely the wrong reason. Assert presence
        separately, where a regression shows up as a plain failure.
        """
        for stid in (ISGF3_CYTOSOL, ISGF3_NUCLEUS, EXPRESSION):
            assert uuids.get(stid), f"{stid} has no uuid in this network"

    @staticmethod
    def _reaches(adj, start: set[str], goal: set[str]) -> bool:
        seen, queue = set(start), list(start)
        while queue:
            for nxt in adj.get(queue.pop(), ()):
                if nxt in goal:
                    return True
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        return False

    def test_nuclear_isgf3_reaches_the_readout(self, forward, uuids):
        """Passes today — the nuclear half of the chain is intact."""
        assert self._reaches(forward, uuids[ISGF3_NUCLEUS], uuids[EXPRESSION]), (
            "nuclear ISGF3 no longer reaches Expression of IFN-induced genes; "
            "the break has moved downstream of nuclear import."
        )

    @pytest.mark.xfail(
        strict=True,
        reason="The cytosolic half is severed at the boundary bridge, so no "
               "perturbation upstream of ISGF3 reaches the readout. This is the "
               "end-to-end symptom of the Tier 2 and Tier 3 failures.",
    )
    def test_cytosolic_isgf3_reaches_the_readout(self, forward, uuids):
        assert self._reaches(forward, uuids[ISGF3_CYTOSOL], uuids[EXPRESSION])
