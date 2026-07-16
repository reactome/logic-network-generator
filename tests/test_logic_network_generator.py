"""Tests for logic_network_generator module."""

from typing import Dict, List, Any
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd

# Add project root to Python path dynamically
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import (
        _assign_uuids,
        _build_uid_index,
        _build_entity_producer_count,
        _canonicalize_registry,
        _emit_boundary_decomposition_edges,
        _register_entity_uuid,
        _get_or_create_entity_uuid,
        _resolve_to_terminal_reactome_ids,
        _resolve_vr_entities,
        export_entity_reaction_proxy_mapping,
    )


class Test_assign_uuids:
    """Tests for _assign_uuids function (position-aware version with union-find)."""

    def test_assigns_new_uuid_for_new_reactome_id(self):
        """Should create a new UUID for a reactome ID not in the registry."""
        entity_uuid_registry: Dict[tuple, str] = {}
        reactome_ids = ["12345"]
        source_reaction_uuid = "source-rxn-uuid"
        target_reaction_uuid = "target-rxn-uuid"

        result = _assign_uuids(reactome_ids, source_reaction_uuid, target_reaction_uuid, entity_uuid_registry)

        assert len(result) == 1
        # Should create entries in registry for both input and output positions
        target_key = ("12345", target_reaction_uuid, "input")
        source_key = ("12345", source_reaction_uuid, "output")
        assert target_key in entity_uuid_registry
        assert source_key in entity_uuid_registry
        # Both should map to same UUID (union-find merged them)
        assert entity_uuid_registry[target_key] == entity_uuid_registry[source_key]
        assert result[0] == entity_uuid_registry[target_key]

    def test_reuses_existing_uuid_for_known_reactome_id_at_same_position(self):
        """Should reuse existing UUID for same reactome ID at same position."""
        existing_uuid = "test-uuid-123"
        source_reaction_uuid = "source-rxn-uuid"
        target_reaction_uuid = "target-rxn-uuid"
        entity_uuid_registry = {
            ("12345", target_reaction_uuid, "input"): existing_uuid,
            ("12345", source_reaction_uuid, "output"): existing_uuid,
        }
        reactome_ids = ["12345"]

        result = _assign_uuids(reactome_ids, source_reaction_uuid, target_reaction_uuid, entity_uuid_registry)

        assert len(result) == 1
        assert result[0] == existing_uuid

    def test_handles_multiple_reactome_ids(self):
        """Should handle multiple reactome IDs correctly at same position."""
        source_reaction_uuid = "source-rxn-uuid"
        target_reaction_uuid = "target-rxn-uuid"
        existing_uuid = "existing-uuid"
        entity_uuid_registry: Dict[tuple, str] = {
            ("12345", target_reaction_uuid, "input"): existing_uuid,
            ("12345", source_reaction_uuid, "output"): existing_uuid,
        }
        reactome_ids = ["12345", "67890", "11111"]

        result = _assign_uuids(reactome_ids, source_reaction_uuid, target_reaction_uuid, entity_uuid_registry)

        assert len(result) == 3
        assert result[0] == existing_uuid  # Reused
        assert result[1] != result[2]  # New UUIDs are different
        assert result[1] != result[0]  # New UUIDs different from existing

    def test_different_positions_get_different_uuids(self):
        """Same reactome ID at different positions should get different UUIDs."""
        entity_uuid_registry: Dict[tuple, str] = {}
        reactome_id = "12345"

        # First position (between reaction1 and reaction2)
        result1 = _assign_uuids([reactome_id], "reaction1-uuid", "reaction2-uuid", entity_uuid_registry)

        # Second position (between reaction3 and reaction4)
        result2 = _assign_uuids([reactome_id], "reaction3-uuid", "reaction4-uuid", entity_uuid_registry)

        # Should have different UUIDs (completely different positions)
        assert result1[0] != result2[0], "Same entity at different positions should have different UUIDs"

    def test_union_find_respects_input_output_roles(self):
        """Entity as input vs output of same reaction should get different UUIDs."""
        entity_uuid_registry: Dict[tuple, str] = {}
        reactome_id = "12345"

        # First edge: reaction1 -> entity -> reaction2 (entity is INPUT to reaction2)
        result1 = _assign_uuids([reactome_id], "reaction1-uuid", "reaction2-uuid", entity_uuid_registry)
        uuid1 = result1[0]

        # Second edge: reaction2 -> entity -> reaction3 (entity is OUTPUT of reaction2)
        result2 = _assign_uuids([reactome_id], "reaction2-uuid", "reaction3-uuid", entity_uuid_registry)
        uuid2 = result2[0]

        # Different roles at same reaction = different positions = different UUIDs
        assert uuid1 != uuid2, "Entity as input vs output of same reaction should have different UUIDs"


class TestEntityProducerCount:
    """Tests for _build_entity_producer_count helper."""

    def test_entity_produced_by_multiple_vrs(self):
        """Entity in output_ids of 2 VRs should have count=2."""
        vr_entities = {
            "vr1": (["A"], ["C", "D"]),
            "vr2": (["B"], ["C", "E"]),
        }
        count = _build_entity_producer_count(vr_entities)
        assert count["C"] == 2
        assert count["D"] == 1
        assert count["E"] == 1

    def test_entity_only_input_not_counted(self):
        """Entity only in input_ids should not appear in count."""
        vr_entities = {
            "vr1": (["A", "B"], ["C"]),
        }
        count = _build_entity_producer_count(vr_entities)
        assert "A" not in count
        assert "B" not in count
        assert count["C"] == 1

    def test_single_producer_returns_one(self):
        """Entity in output_ids of 1 VR should have count=1."""
        vr_entities = {
            "vr1": (["A"], ["X"]),
            "vr2": (["B"], ["Y"]),
        }
        count = _build_entity_producer_count(vr_entities)
        assert count["X"] == 1
        assert count["Y"] == 1


class TestUidIndex:
    """Tests for cached decomposition UID resolution."""

    def test_stringified_nulls_are_not_treated_as_entities(self):
        """Cached CSV reloads can contain literal 'None' strings; skip them."""
        decomposed = pd.DataFrame(
            [
                {
                    "uid": "hash-1",
                    "input_or_output_uid": "None",
                    "input_or_output_reactome_id": "R-HSA-1",
                    "stoichiometry": 2,
                },
                {
                    "uid": "hash-1",
                    "input_or_output_uid": "",
                    "input_or_output_reactome_id": "nan",
                    "stoichiometry": 1,
                },
                {
                    "uid": "hash-1",
                    "input_or_output_uid": "hash-2",
                    "input_or_output_reactome_id": "<NA>",
                    "stoichiometry": 3,
                },
                {
                    "uid": "hash-2",
                    "input_or_output_uid": "NaN",
                    "input_or_output_reactome_id": "R-HSA-2",
                    "stoichiometry": 4,
                },
            ]
        )

        uid_index = _build_uid_index(decomposed)

        assert uid_index["hash-1"][0] == ["hash-2"]
        assert uid_index["hash-1"][1] == ["R-HSA-1"]
        assert "None" not in uid_index["hash-1"][2]
        assert "nan" not in uid_index["hash-1"][2]

        resolved = _resolve_to_terminal_reactome_ids(uid_index, "hash-1")

        assert resolved == {"R-HSA-1": 2, "R-HSA-2": 12}


class TestInterReactionConnectivity:
    """Tests for inter-reaction entity UUID connectivity (3-phase approach).

    Verifies that entities shared between reactions get merged UUIDs,
    while disconnected entities remain separate.
    """

    def test_two_reactions_share_entity_uuid(self):
        """Entity shared as output of VR1 and input of VR2 should get one UUID."""
        registry: Dict[tuple, str] = {}
        unions: Dict[str, str] = {}

        # Phase 1: Register
        _register_entity_uuid("A", "vr1", "output", registry)
        _register_entity_uuid("A", "vr2", "input", registry)

        # Should start as different UUIDs
        assert registry[("A", "vr1", "output")] != registry[("A", "vr2", "input")]

        # Phase 2: Merge
        _get_or_create_entity_uuid("A", "vr1", "vr2", registry, unions)
        _canonicalize_registry(registry, unions)

        # Should now share the same UUID
        assert registry[("A", "vr1", "output")] == registry[("A", "vr2", "input")]

    def test_three_reaction_chain(self):
        """VR1→A→VR2→B→VR3: A and B have separate merged UUIDs."""
        registry: Dict[tuple, str] = {}
        unions: Dict[str, str] = {}

        # Phase 1: Register all entities
        _register_entity_uuid("A", "vr1", "output", registry)
        _register_entity_uuid("A", "vr2", "input", registry)
        _register_entity_uuid("B", "vr2", "output", registry)
        _register_entity_uuid("B", "vr3", "input", registry)

        # Phase 2: Merge connections
        _get_or_create_entity_uuid("A", "vr1", "vr2", registry, unions)
        _get_or_create_entity_uuid("B", "vr2", "vr3", registry, unions)
        _canonicalize_registry(registry, unions)

        uuid_a = registry[("A", "vr1", "output")]
        uuid_b = registry[("B", "vr2", "output")]

        # A and B should have different UUIDs
        assert uuid_a != uuid_b

        # A consistent across VR1 output and VR2 input
        assert registry[("A", "vr1", "output")] == registry[("A", "vr2", "input")]

        # B consistent across VR2 output and VR3 input
        assert registry[("B", "vr2", "output")] == registry[("B", "vr3", "input")]

    def test_no_spurious_keys(self):
        """_register_entity_uuid should create only one key per call."""
        registry: Dict[tuple, str] = {}

        _register_entity_uuid("A", "vr1", "input", registry)

        assert len(registry) == 1
        assert ("A", "vr1", "input") in registry
        assert ("A", "vr1", "output") not in registry

    def test_disconnected_reactions_different_uuids(self):
        """Same entity in unconnected reactions should have different UUIDs."""
        registry: Dict[tuple, str] = {}

        _register_entity_uuid("A", "vr1", "output", registry)
        _register_entity_uuid("A", "vr3", "input", registry)

        # No Phase 2 merge — they're disconnected
        assert registry[("A", "vr1", "output")] != registry[("A", "vr3", "input")]

    def test_multi_source_convergence(self):
        """VR1→A→VR2 and VR3→A→VR2 should all merge to same UUID."""
        registry: Dict[tuple, str] = {}
        unions: Dict[str, str] = {}

        # Phase 1: Register
        _register_entity_uuid("A", "vr1", "output", registry)
        _register_entity_uuid("A", "vr3", "output", registry)
        _register_entity_uuid("A", "vr2", "input", registry)

        # Phase 2: Both VR1 and VR3 feed A into VR2
        _get_or_create_entity_uuid("A", "vr1", "vr2", registry, unions)
        _get_or_create_entity_uuid("A", "vr3", "vr2", registry, unions)
        _canonicalize_registry(registry, unions)

        uuid_from_vr1 = registry[("A", "vr1", "output")]
        uuid_from_vr3 = registry[("A", "vr3", "output")]
        uuid_at_vr2 = registry[("A", "vr2", "input")]

        # All three should share the same UUID
        assert uuid_from_vr1 == uuid_at_vr2
        assert uuid_from_vr3 == uuid_at_vr2

    def test_no_duplicate_edges(self, monkeypatch):
        """_resolve_vr_entities must emit each node once (Set-deduped).

        Since set-variant emission, node identities come from the reaction's
        annotated input/output entities (mapped to nodes), and duplicates are
        collapsed via a set. We mock the two Neo4j calls the resolver makes:
        the reaction's annotated entities and their labels (simple proteins).
        """
        import src.logic_network_generator as m
        from src import neo4j_connector

        # Reaction "1" annotates one simple-protein input and one output.
        monkeypatch.setattr(
            neo4j_connector, "get_reaction_input_output_ids",
            lambda rid, io: {"9933417"} if io == "input" else {"12345"},
        )
        # Both entities are simple (not complexes/sets) -> mapped to themselves.
        monkeypatch.setattr(
            neo4j_connector, "get_labels",
            lambda e: ["EntityWithAccessionedSequence"],
        )

        # input_hash resolves to member "9933417" via two convergent paths.
        uid_index = {
            "vr1-input": (["nested-1", "nested-2"], set(), {}),
            "nested-1": ([], {"9933417"}, {"9933417": 1}),
            "nested-2": ([], {"9933417"}, {"9933417": 1}),
            "vr1-output": ([], {"12345"}, {"12345": 1}),
        }
        reaction_id_map = pd.DataFrame({
            "uid": ["vr1"],
            "input_hash": ["vr1-input"],
            "output_hash": ["vr1-output"],
            "reactome_id": [1],
        })

        vr_entities = _resolve_vr_entities(reaction_id_map, uid_index)
        input_ids, output_ids, _in_stoich, _out_stoich = vr_entities["vr1"]

        assert input_ids == ["9933417"], f"expected one deduped input, got {input_ids}"
        assert output_ids == ["12345"], f"expected one output, got {output_ids}"

    def test_root_input_same_entity_gets_one_uuid(self):
        """Root input entity appearing at multiple reactions should share one UUID."""
        registry: Dict[tuple, str] = {}
        root_input_eids = {"A"}
        root_input_cache: Dict[str, str] = {}

        _register_entity_uuid("A", "vr1", "input", registry,
                              root_input_eids, root_input_cache)
        _register_entity_uuid("A", "vr3", "input", registry,
                              root_input_eids, root_input_cache)

        assert registry[("A", "vr1", "input")] == registry[("A", "vr3", "input")]

    def test_terminal_output_same_entity_gets_one_uuid(self):
        """Terminal output entity appearing at multiple reactions should share one UUID."""
        registry: Dict[tuple, str] = {}
        terminal_output_eids = {"B"}
        terminal_output_cache: Dict[str, str] = {}

        _register_entity_uuid("B", "vr1", "output", registry,
                              terminal_output_eids, terminal_output_cache)
        _register_entity_uuid("B", "vr2", "output", registry,
                              terminal_output_eids, terminal_output_cache)

        assert registry[("B", "vr1", "output")] == registry[("B", "vr2", "output")]

    def test_root_and_terminal_same_entity_different_uuids(self):
        """Entity that is both root input and terminal output should get separate UUIDs."""
        registry: Dict[tuple, str] = {}
        root_input_eids = {"A"}
        terminal_output_eids = {"A"}
        root_cache: Dict[str, str] = {}
        terminal_cache: Dict[str, str] = {}

        _register_entity_uuid("A", "vr1", "input", registry,
                              root_input_eids, root_cache)
        _register_entity_uuid("A", "vr2", "output", registry,
                              terminal_output_eids, terminal_cache)

        # Different caches → different UUIDs
        assert registry[("A", "vr1", "input")] != registry[("A", "vr2", "output")]

    def test_non_boundary_entity_gets_separate_uuids(self):
        """Entity not in boundary sets should get normal per-position UUIDs."""
        registry: Dict[tuple, str] = {}
        root_input_eids = {"X"}  # "A" is NOT a boundary entity
        root_cache: Dict[str, str] = {}

        _register_entity_uuid("A", "vr1", "input", registry,
                              root_input_eids, root_cache)
        _register_entity_uuid("A", "vr2", "input", registry,
                              root_input_eids, root_cache)

        # "A" is not in root_input_eids, so it gets separate UUIDs
        assert registry[("A", "vr1", "input")] != registry[("A", "vr2", "input")]


class TestBoundaryLeavesReuseExistingUUIDs:
    """Boundary expansion must NOT mint a fresh UUID for a leaf if that
    leaf's stId already has a UUID elsewhere in the network. Otherwise
    perturbing 'MDM2 the regulator' wouldn't propagate to 'MDM2 the
    boundary leaf' — they'd be disconnected duplicate nodes for the
    same biological entity, which kills the perturbation use case
    boundary expansion exists to enable.
    """

    def _root_edge(self, complex_uuid):
        # Seed an edge that makes complex_uuid a root input: a source that is
        # never a target (its target is an unmapped reaction node).
        return [{"source_id": complex_uuid, "target_id": "u-reaction",
                 "pos_neg": "pos", "and_or": "and",
                 "edge_type": "input", "stoichiometry": 1}]

    def test_leaf_reuses_uuid_when_entity_already_in_registry(self):
        """If MDM2 already has UUID U_existing in reactome_id_to_uuid (e.g.
        because it's a regular VR input or a regulator elsewhere), the
        boundary expansion of root-input MDM2:TP53 must use U_existing for
        the MDM2 leaf — not a fresh one.
        """
        existing_mdm2_uuid = "u-existing-mdm2"
        complex_uuid = "u-complex"
        reactome_id_to_uuid: Dict[str, str] = {
            existing_mdm2_uuid: "MDM2",  # MDM2 already has a UUID elsewhere
            complex_uuid: "MDM2:TP53",
        }
        edges: List[Dict[str, Any]] = self._root_edge(complex_uuid)

        with patch('src.neo4j_connector.get_labels',
                   return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"MDM2", "TP53"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                reactome_id_to_uuid=reactome_id_to_uuid,
            )

        mdm2_assembly = [
            e for e in edges
            if e["edge_type"] == "assembly"
            and e["target_id"] == complex_uuid
            and reactome_id_to_uuid.get(e["source_id"]) == "MDM2"
        ]
        assert len(mdm2_assembly) == 1
        assert mdm2_assembly[0]["source_id"] == existing_mdm2_uuid, (
            "Boundary leaf must reuse the existing UUID for its stId, "
            "not mint a fresh one. See docs/DESIGN_DECISIONS.md."
        )

    def test_leaf_mints_fresh_uuid_when_entity_is_new_to_network(self):
        """When the leaf's stId is not in reactome_id_to_uuid, a fresh
        UUID is minted and recorded so future passes can find it.
        """
        complex_uuid = "u-complex"
        reactome_id_to_uuid: Dict[str, str] = {complex_uuid: "C"}
        edges: List[Dict[str, Any]] = self._root_edge(complex_uuid)

        with patch('src.neo4j_connector.get_labels',
                   return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"L1", "L2"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                reactome_id_to_uuid=reactome_id_to_uuid,
            )

        new_uuids = [u for u, sid in reactome_id_to_uuid.items() if sid in {"L1", "L2"}]
        assert len(new_uuids) == 2
        for u in new_uuids:
            assert any(e["source_id"] == u and e["edge_type"] == "assembly" for e in edges)

    def test_same_complex_at_many_roots_shares_one_member_node(self):
        """A complex appearing as a root input at N positions yields ONE member
        node with N assembly edges — members are not duplicated (the user's
        requirement)."""
        c1, c2, c3 = "u-c1", "u-c2", "u-c3"   # three root-input occurrences
        reactome_id_to_uuid: Dict[str, str] = {
            c1: "MDC1c", c2: "MDC1c", c3: "MDC1c",  # same complex stId, 3 UUIDs
        }
        edges: List[Dict[str, Any]] = []
        for c in (c1, c2, c3):
            edges += self._root_edge(c)

        with patch('src.neo4j_connector.get_labels', return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"MDC1"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                reactome_id_to_uuid=reactome_id_to_uuid,
            )

        asm = [e for e in edges if e["edge_type"] == "assembly"]
        member_uuids = {e["source_id"] for e in asm}
        complex_targets = {e["target_id"] for e in asm}
        assert len(member_uuids) == 1, "MDC1 member must be a single shared node"
        assert complex_targets == {c1, c2, c3}, "one assembly edge to each occurrence"

    def test_intermediate_occurrence_not_decomposed(self):
        """A complex that is produced and consumed (intermediate) is left intact."""
        cx = "u-cx"
        reactome_id_to_uuid: Dict[str, str] = {cx: "INT"}
        # cx is both a target (produced) and a source (consumed) → intermediate.
        edges: List[Dict[str, Any]] = [
            {"source_id": "u-rxn1", "target_id": cx, "pos_neg": "pos",
             "and_or": "and", "edge_type": "output", "stoichiometry": 1},
            {"source_id": cx, "target_id": "u-rxn2", "pos_neg": "pos",
             "and_or": "and", "edge_type": "input", "stoichiometry": 1},
        ]
        with patch('src.neo4j_connector.get_labels', return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"M1", "M2"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                reactome_id_to_uuid=reactome_id_to_uuid,
            )
        assert not any(e["edge_type"] in ("assembly", "dissociation") for e in edges), \
            "intermediate complex must not be decomposed"

    def test_dissociation_members_are_separate_readout_sinks(self):
        """Terminal-output complex members come out as FRESH sink nodes — not the
        member's functional node — with no outgoing edges, so the complex can't
        inject its activity into the member's other roles."""
        functional_mdc1 = "u-functional-mdc1"
        term_complex = "u-termc"
        reactome_id_to_uuid: Dict[str, str] = {
            functional_mdc1: "MDC1",        # MDC1 already a functional node
            term_complex: "MDC1:partner",   # the terminal-output complex
        }
        edges: List[Dict[str, Any]] = [
            # functional MDC1 is intermediate (produced + consumed) → left intact
            {"source_id": "u-rxn0", "target_id": functional_mdc1, "pos_neg": "pos",
             "and_or": "and", "edge_type": "output", "stoichiometry": 1},
            {"source_id": functional_mdc1, "target_id": "u-rxn2", "pos_neg": "pos",
             "and_or": "and", "edge_type": "input", "stoichiometry": 1},
            # term_complex is produced, never consumed → terminal output
            {"source_id": "u-rxn", "target_id": term_complex, "pos_neg": "pos",
             "and_or": "and", "edge_type": "output", "stoichiometry": 1},
        ]
        with patch('src.neo4j_connector.get_labels', return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"MDC1", "PARTNER"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                reactome_id_to_uuid=reactome_id_to_uuid,
            )
        diss = [e for e in edges if e["edge_type"] == "dissociation"]
        assert len(diss) == 2, "one readout sink per member"
        sinks = {e["target_id"] for e in diss}
        assert functional_mdc1 not in sinks, \
            "dissociation must NOT reuse the member's functional node"
        all_sources = {e["source_id"] for e in edges}
        for s in sinks:
            assert s not in all_sources, "readout sink must have no outgoing edges"
            assert reactome_id_to_uuid[s] in {"MDC1", "PARTNER"}, \
                "sink must carry the member's stId for measurement"


class TestEntityReactionProxyMapping:
    """Tests for export_entity_reaction_proxy_mapping.

    A curated species (often a Complex containing an EntitySet) gets expanded
    into virtual variants during generation, so its own stId never appears in
    stid_to_uuid_mapping.csv. This export restores addressability by pointing
    the species at the UUIDs of the reaction that produces it.
    """

    def _write(self, tmp_path, network, reaction_id_map, reactome_id_to_uuid,
               participating, entity_reactions):
        out = tmp_path / "proxy.csv"
        with patch('src.neo4j_connector.get_pathway_participating_entities',
                   return_value=participating), \
             patch('src.neo4j_connector.get_pathway_entity_reactions',
                   return_value=entity_reactions):
            export_entity_reaction_proxy_mapping(
                network, reaction_id_map, reactome_id_to_uuid,
                "R-HSA-1", str(out),
            )
        return pd.read_csv(out)

    # Fake UUIDs must contain a dash: export_*_mapping detects the dict
    # orientation with the same `'-' in key` heuristic the rest of the module
    # uses, and real position-aware UUIDs always have dashes.
    def test_missing_complex_maps_to_producing_reaction(self, tmp_path):
        # Network: reaction R1 (uuid u-r1) outputs the expanded variant nodes;
        # the parent complex C is absent from the mapping.
        network = pd.DataFrame({
            "source_id": ["u-in", "u-r1"],
            "target_id": ["u-r1", "u-out"],
        })
        reaction_id_map = pd.DataFrame({"uid": ["u-r1"], "reactome_id": ["R-HSA-R1"]})
        # Mapping has the reaction and some leaves, but NOT complex C.
        reactome_id_to_uuid = {"u-r1": "R-HSA-R1", "u-in": "R-HSA-IN", "u-out": "R-HSA-OUT"}

        df = self._write(
            tmp_path, network, reaction_id_map, reactome_id_to_uuid,
            participating={"R-HSA-C", "R-HSA-IN", "R-HSA-OUT"},
            entity_reactions={"R-HSA-C": {"output": ["R-HSA-R1"]}},
        )
        rows = df[df["entity_stable_id"] == "R-HSA-C"]
        assert list(rows["proxy_uuid"]) == ["u-r1"]
        assert set(rows["proxy_role"]) == {"producing"}

    def test_entity_already_in_mapping_is_skipped(self, tmp_path):
        network = pd.DataFrame({"source_id": ["u-r1"], "target_id": ["u-out"]})
        reaction_id_map = pd.DataFrame({"uid": ["u-r1"], "reactome_id": ["R-HSA-R1"]})
        # R-HSA-OUT is directly present, so it must not be proxied.
        reactome_id_to_uuid = {"u-r1": "R-HSA-R1", "u-out": "R-HSA-OUT"}
        df = self._write(
            tmp_path, network, reaction_id_map, reactome_id_to_uuid,
            participating={"R-HSA-OUT"},
            entity_reactions={"R-HSA-OUT": {"output": ["R-HSA-R1"]}},
        )
        assert df.empty or "R-HSA-OUT" not in set(df["entity_stable_id"])

    def test_consuming_fallback_when_no_producer(self, tmp_path):
        network = pd.DataFrame({"source_id": ["u-r1"], "target_id": ["u-out"]})
        reaction_id_map = pd.DataFrame({"uid": ["u-r1"], "reactome_id": ["R-HSA-R1"]})
        reactome_id_to_uuid = {"u-r1": "R-HSA-R1", "u-out": "R-HSA-OUT"}
        # C is only ever consumed (input) — no producing reaction in-pathway.
        df = self._write(
            tmp_path, network, reaction_id_map, reactome_id_to_uuid,
            participating={"R-HSA-C"},
            entity_reactions={"R-HSA-C": {"input": ["R-HSA-R1"]}},
        )
        rows = df[df["entity_stable_id"] == "R-HSA-C"]
        assert list(rows["proxy_uuid"]) == ["u-r1"]
        assert set(rows["proxy_role"]) == {"consuming"}

    def test_producing_preferred_over_consuming(self, tmp_path):
        network = pd.DataFrame({
            "source_id": ["u-r1", "u-r2"], "target_id": ["u-x", "u-y"],
        })
        reaction_id_map = pd.DataFrame(
            {"uid": ["u-r1", "u-r2"], "reactome_id": ["R-HSA-R1", "R-HSA-R2"]})
        reactome_id_to_uuid = {"u-r1": "R-HSA-R1", "u-r2": "R-HSA-R2"}
        # C is produced by R1 and consumed by R2 — only the producer should win.
        df = self._write(
            tmp_path, network, reaction_id_map, reactome_id_to_uuid,
            participating={"R-HSA-C"},
            entity_reactions={"R-HSA-C": {"output": ["R-HSA-R1"], "input": ["R-HSA-R2"]}},
        )
        rows = df[df["entity_stable_id"] == "R-HSA-C"]
        assert list(rows["proxy_uuid"]) == ["u-r1"]
        assert set(rows["proxy_role"]) == {"producing"}

    def test_reaction_not_in_network_yields_no_row(self, tmp_path):
        # The producing reaction exists in Reactome but its UUID never made it
        # into the network (e.g. dropped for having no I/O) — nothing to proxy.
        network = pd.DataFrame({"source_id": ["u-a"], "target_id": ["u-b"]})
        reaction_id_map = pd.DataFrame({"uid": [], "reactome_id": []})
        reactome_id_to_uuid = {"u-a": "R-HSA-A", "u-b": "R-HSA-B"}
        df = self._write(
            tmp_path, network, reaction_id_map, reactome_id_to_uuid,
            participating={"R-HSA-C"},
            entity_reactions={"R-HSA-C": {"output": ["R-HSA-MISSING"]}},
        )
        assert df.empty or "R-HSA-C" not in set(df["entity_stable_id"])


class TestRegulatorVariantCap:
    """_expand_complex_variants (the negative-regulator decomposition path)
    must bundle a complex into ONE opaque node when its set-variant cartesian
    exceeds MAX_VARIANTS (issue #40). Without this a single inhibitor complex
    with large internal EntitySets fans out into ~144k variant nodes, each
    wired to every reaction it inhibits — RAF/MAP kinase hit 4.76M edges.
    """

    def _wire_mocks(self, monkeypatch):
        import src.logic_network_generator as m
        from src import neo4j_connector as nc
        # Complex C = two EntitySets of 3 members each -> 3x3 = 9 variants.
        monkeypatch.setattr(nc, "get_labels", lambda e: {
            "C": ["Complex"], "S1": ["EntitySet"], "S2": ["EntitySet"],
        }.get(e, ["EntityWithAccessionedSequence"]))
        monkeypatch.setattr(nc, "get_complex_components",
                            lambda e: {"S1": 1, "S2": 1} if e == "C" else {})
        monkeypatch.setattr(nc, "get_set_members", lambda e: {
            "S1": ["A", "B", "D"], "S2": ["E", "F", "G"]}.get(e, []))
        monkeypatch.setattr(m, "_complex_contains_entity_set", lambda e: e == "C")
        return m

    def test_over_cap_bundles_to_single_node(self, monkeypatch):
        m = self._wire_mocks(monkeypatch)
        monkeypatch.setattr(m, "MAX_VARIANTS", 4)
        result = m._expand_complex_variants("C")
        assert result == [("C", 1)], (
            f"9 variants > cap 4 must bundle to one node, got {result}"
        )

    def test_under_cap_still_expands(self, monkeypatch):
        m = self._wire_mocks(monkeypatch)
        monkeypatch.setattr(m, "MAX_VARIANTS", 100)
        result = m._expand_complex_variants("C")
        assert len(result) == 9, (
            f"3x3 cartesian under cap should yield 9 variants, got {len(result)}"
        )
        assert all(vid.startswith("C::variant::") for vid, _ in result)


class TestMatchingLeavesGranularity:
    """_matching_leaves mirrors the MATCHING layer's decomposition depth so a
    set participant lines up with the terminal reactome-ids a VR resolved to.

    The bug it fixes: get_terminal_components collapses modified species down to
    their base reference protein (p-ERK -> ERK), while the matching layer keeps
    the modified form. Intersecting the two came up empty, so EntitySets fell
    back to bare, producer-less nodes. _matching_leaves keeps Complexes atomic
    (unless they contain a set) — the same rule break_apart_entity uses.
    """

    def _clear(self):
        import src.logic_network_generator as m
        m._matching_leaves_cache.clear()

    def test_set_of_atomic_complexes_stops_at_the_complexes(self, monkeypatch):
        import src.logic_network_generator as m
        from src import neo4j_connector as nc
        self._clear()
        # Set S -> {C1, C2}; C1,C2 are Complexes WITHOUT internal sets (atomic).
        monkeypatch.setattr(nc, "get_labels", lambda e: {
            "S": ["EntitySet", "DefinedSet"], "C1": ["Complex"], "C2": ["Complex"],
        }.get(e, ["EntityWithAccessionedSequence"]))
        monkeypatch.setattr(nc, "get_set_members",
                            lambda e: {"C1", "C2"} if e == "S" else set())
        monkeypatch.setattr(nc, "get_complex_components",
                            lambda e: {"P1": 1, "P2": 1} if e in ("C1", "C2") else {})
        monkeypatch.setattr(m, "_complex_contains_entity_set", lambda e: False)
        result = m._matching_leaves("S")
        assert result == frozenset({"C1", "C2"}), (
            "a set of atomic complexes must resolve to the complexes (matching "
            f"granularity), not their proteins; got {set(result)}"
        )

    def test_complex_with_internal_set_decomposes(self, monkeypatch):
        import src.logic_network_generator as m
        from src import neo4j_connector as nc
        self._clear()
        # Complex C contains set S={A,B} plus protein P -> decompose to {A,B,P}.
        monkeypatch.setattr(nc, "get_labels", lambda e: {
            "C": ["Complex"], "S": ["EntitySet"],
        }.get(e, ["EntityWithAccessionedSequence"]))
        monkeypatch.setattr(nc, "get_complex_components",
                            lambda e: {"S": 1, "P": 1} if e == "C" else {})
        monkeypatch.setattr(nc, "get_set_members",
                            lambda e: {"A", "B"} if e == "S" else set())
        monkeypatch.setattr(m, "_complex_contains_entity_set", lambda e: e == "C")
        result = m._matching_leaves("C")
        assert result == frozenset({"A", "B", "P"}), (
            f"complex-with-set must decompose to members, got {set(result)}"
        )


class TestComplexAsNode:
    """LNG_COMPLEX_AS_NODE (default on): a Complex is ONE node = its plain stId
    even when it contains an internal EntitySet, so the produced-complex node id
    matches the catalyst/regulator node id and they unify (fixes the UUID silo
    that severed produce->catalyze chains, e.g. ATM->p53:FAS-gene complex->FAS)."""

    def test_complex_with_set_emits_plain_stid(self, monkeypatch):
        import src.logic_network_generator as m
        from src import neo4j_connector as nc
        monkeypatch.setattr(nc, "get_labels", lambda e: ["Complex"] if e == "C" else ["EntitySet"])
        monkeypatch.setattr(m, "_complex_contains_entity_set", lambda e: True)
        # default-on: complex-with-set -> plain stId, never a ::variant:: node
        result = m._map_annotated_entity_to_nodes("C", {"C", "A", "B"})
        assert result == {"C"}, f"complex-as-node must emit plain stId, got {result}"

    def test_flag_off_restores_variant_behavior(self, monkeypatch):
        import src.logic_network_generator as m
        from src import neo4j_connector as nc
        monkeypatch.setenv("LNG_COMPLEX_AS_NODE", "0")
        monkeypatch.setattr(nc, "get_labels", lambda e: ["Complex"])
        monkeypatch.setattr(m, "_complex_contains_entity_set", lambda e: True)
        monkeypatch.setattr(m, "_complex_variant_leafsets", lambda e: [frozenset({"A", "B"})])
        m._variant_capped.discard("C")
        result = m._map_annotated_entity_to_nodes("C", {"A", "B"})
        assert any("::variant::" in r for r in result), f"flag off must variant-split, got {result}"
