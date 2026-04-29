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
        _build_entity_producer_count,
        _canonicalize_registry,
        _emit_boundary_decomposition_edges,
        _register_entity_uuid,
        _get_or_create_entity_uuid,
        _resolve_vr_entities,
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

    def test_no_duplicate_edges(self):
        """Duplicate terminal IDs from decomposition should not produce duplicate edges.

        When multiple decomposition paths converge on the same terminal Reactome ID,
        _resolve_to_terminal_reactome_ids returns duplicates. _resolve_vr_entities
        must deduplicate them so Phase 3 doesn't create duplicate edges.
        """
        # Build a uid_index where hash "vr1-input" resolves to terminal ID "9933417"
        # via two different nested paths, producing duplicates without dedup.
        # uid_index maps hash -> (nested_uids, terminal_ids, stoich_map)
        uid_index = {
            "vr1-input": (["nested-1", "nested-2"], set(), {}),  # two nested paths, no direct terminals
            "nested-1": ([], {"9933417"}, {"9933417": 1}),  # both nested paths resolve to same terminal
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

        input_ids, output_ids, input_stoich, output_stoich = vr_entities["vr1"]

        # _resolve_to_terminal_reactome_ids now returns dict (deduped by key),
        # but stoichiometry accumulates: 1 + 1 = 2 from two nested paths
        assert len(input_ids) == 1, (
            f"Expected 1 unique input ID, got {len(input_ids)}: {input_ids}"
        )
        assert input_ids[0] == "9933417"
        assert input_stoich["9933417"] == 2  # stoichiometry adds: 1 from nested-1 + 1 from nested-2
        assert len(output_ids) == 1

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

    def test_leaf_reuses_uuid_when_entity_already_in_registry(self):
        """If MDM2 already has UUID U_existing in reactome_id_to_uuid (e.g.
        because it's a regular VR input or a regulator elsewhere), the
        boundary expansion of MDM2:TP53 must use U_existing for the MDM2
        leaf — not a fresh one.
        """
        existing_mdm2_uuid = "u-existing-mdm2"
        complex_uuid = "u-complex"
        reactome_id_to_uuid: Dict[str, str] = {
            existing_mdm2_uuid: "MDM2",  # MDM2 already has a UUID elsewhere
            complex_uuid: "MDM2:TP53",
        }
        edges: List[Dict[str, Any]] = []

        with patch('src.neo4j_connector.get_labels',
                   return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"MDM2", "TP53"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                root_input_eids={"MDM2:TP53"},
                terminal_output_eids=set(),
                root_input_uuid_cache={"MDM2:TP53": complex_uuid},
                terminal_output_uuid_cache={},
                reactome_id_to_uuid=reactome_id_to_uuid,
            )

        # Find the assembly edge whose target is the complex and whose
        # source maps back to MDM2.
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
        edges: List[Dict[str, Any]] = []

        with patch('src.neo4j_connector.get_labels',
                   return_value=["Complex"]), \
             patch('src.logic_network_generator.get_terminal_components',
                   return_value={"L1", "L2"}):
            _emit_boundary_decomposition_edges(
                pathway_logic_network_data=edges,
                root_input_eids={"C"},
                terminal_output_eids=set(),
                root_input_uuid_cache={"C": complex_uuid},
                terminal_output_uuid_cache={},
                reactome_id_to_uuid=reactome_id_to_uuid,
            )

        # Two new leaf UUIDs added to the mapping
        new_uuids = [u for u, sid in reactome_id_to_uuid.items() if sid in {"L1", "L2"}]
        assert len(new_uuids) == 2
        # Each fresh UUID is also used as a source on an assembly edge
        for u in new_uuids:
            assert any(e["source_id"] == u and e["edge_type"] == "assembly" for e in edges)
