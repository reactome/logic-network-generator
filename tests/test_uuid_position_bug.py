"""Test for UUID position-awareness.

This test verifies that the same Reactome entity appearing at different
positions in a pathway receives different UUIDs in the logic network.

The current implementation uses union-find logic with
(entity_dbId, reaction_uuid, role) tuples as keys to ensure entities
at different pathway positions get different UUIDs.
"""

import uuid
import pytest
import sys
from pathlib import Path
from unittest.mock import patch

# Add project root to Python path dynamically
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import _assign_uuids, _get_or_create_entity_uuid


def test_same_entity_different_reactions_get_different_uuids():
    """Test that the same entity in different reaction contexts gets different UUIDs.

    When entity 179838 is an output of reaction A and input to reaction B,
    it should get a different UUID than when it connects reaction C to reaction D.
    """
    entity_uuid_registry = {}

    # Entity 179838 connecting reaction_A -> reaction_B
    reaction_a_uuid = str(uuid.uuid4())
    reaction_b_uuid = str(uuid.uuid4())

    uuid1 = _get_or_create_entity_uuid(
        179838, reaction_a_uuid, reaction_b_uuid, entity_uuid_registry
    )

    # Same entity 179838 connecting reaction_C -> reaction_D
    reaction_c_uuid = str(uuid.uuid4())
    reaction_d_uuid = str(uuid.uuid4())

    uuid2 = _get_or_create_entity_uuid(
        179838, reaction_c_uuid, reaction_d_uuid, entity_uuid_registry
    )

    # Different reaction contexts should produce different UUIDs
    assert uuid1 != uuid2, (
        f"Entity 179838 in different reaction contexts should have DIFFERENT UUIDs.\n"
        f"Context 1 ({reaction_a_uuid[:8]}... -> {reaction_b_uuid[:8]}...): {uuid1}\n"
        f"Context 2 ({reaction_c_uuid[:8]}... -> {reaction_d_uuid[:8]}...): {uuid2}"
    )


def test_same_entity_same_connection_gets_same_uuid():
    """Test that the same entity in the same reaction context gets the same UUID.

    When entity 179838 connects reaction_A output to reaction_B input,
    calling again with the same context should return the same UUID.
    """
    entity_uuid_registry = {}

    reaction_a_uuid = str(uuid.uuid4())
    reaction_b_uuid = str(uuid.uuid4())

    uuid1 = _get_or_create_entity_uuid(
        179838, reaction_a_uuid, reaction_b_uuid, entity_uuid_registry
    )
    uuid2 = _get_or_create_entity_uuid(
        179838, reaction_a_uuid, reaction_b_uuid, entity_uuid_registry
    )

    assert uuid1 == uuid2, (
        f"Same entity in same context should get the SAME UUID.\n"
        f"First call: {uuid1}\nSecond call: {uuid2}"
    )


def test_entity_different_roles_at_same_reaction_get_different_uuids():
    """Test that entity at different roles (input vs output) of the same reaction gets different UUIDs.

    The current implementation uses (entity_dbId, reaction_uuid, role) tuples.
    Entity 179838 as input to reaction_B (from A->B) has a different position
    than entity 179838 as output of reaction_B (from B->C), so they get
    different UUIDs.
    """
    entity_uuid_registry = {}

    reaction_a_uuid = str(uuid.uuid4())
    reaction_b_uuid = str(uuid.uuid4())
    reaction_c_uuid = str(uuid.uuid4())

    # Entity connects A -> B (entity is input to B)
    uuid_ab = _get_or_create_entity_uuid(
        179838, reaction_a_uuid, reaction_b_uuid, entity_uuid_registry
    )

    # Same entity connects B -> C (entity is output of B)
    uuid_bc = _get_or_create_entity_uuid(
        179838, reaction_b_uuid, reaction_c_uuid, entity_uuid_registry
    )

    # Different roles at reaction_b: "input" vs "output" are different positions
    assert uuid_ab != uuid_bc, (
        f"Entity at different roles of same reaction should have DIFFERENT UUIDs.\n"
        f"A->B (input to B): {uuid_ab}\nB->C (output of B): {uuid_bc}"
    )


def test_assign_uuids_batch():
    """Test _assign_uuids assigns UUIDs to multiple entities in batch."""
    entity_uuid_registry = {}

    source_uuid = str(uuid.uuid4())
    target_uuid = str(uuid.uuid4())

    reactome_ids = [179838, 1002, 54321]

    uuids = _assign_uuids(reactome_ids, source_uuid, target_uuid, entity_uuid_registry)

    assert len(uuids) == 3, "Should assign UUID to each entity"
    assert len(set(uuids)) == 3, "Different entities should get different UUIDs"


def test_different_entities_same_context_get_different_uuids():
    """Test that different entities in the same reaction context get different UUIDs."""
    entity_uuid_registry = {}

    reaction_a_uuid = str(uuid.uuid4())
    reaction_b_uuid = str(uuid.uuid4())

    uuid_entity1 = _get_or_create_entity_uuid(
        179838, reaction_a_uuid, reaction_b_uuid, entity_uuid_registry
    )
    uuid_entity2 = _get_or_create_entity_uuid(
        1002, reaction_a_uuid, reaction_b_uuid, entity_uuid_registry
    )

    assert uuid_entity1 != uuid_entity2, (
        f"Different entities should have different UUIDs even in same context.\n"
        f"Entity 179838: {uuid_entity1}\nEntity 1002: {uuid_entity2}"
    )


def test_full_scenario_entity_at_three_positions():
    """Test entity appearing at 3 independent pathway positions.

    Entity 179838 appears at:
    - Position 1: reaction_A -> reaction_B
    - Position 2: reaction_C -> reaction_D
    - Position 3: reaction_E -> reaction_F

    All three should get DIFFERENT UUIDs since they are at different pathway positions.
    """
    entity_uuid_registry = {}

    # Create 6 unique reactions
    reactions = [str(uuid.uuid4()) for _ in range(6)]

    uuid_pos1 = _get_or_create_entity_uuid(179838, reactions[0], reactions[1], entity_uuid_registry)
    uuid_pos2 = _get_or_create_entity_uuid(179838, reactions[2], reactions[3], entity_uuid_registry)
    uuid_pos3 = _get_or_create_entity_uuid(179838, reactions[4], reactions[5], entity_uuid_registry)

    assert uuid_pos1 != uuid_pos2, "Positions 1 & 2 should have DIFFERENT UUIDs"
    assert uuid_pos1 != uuid_pos3, "Positions 1 & 3 should have DIFFERENT UUIDs"
    assert uuid_pos2 != uuid_pos3, "Positions 2 & 3 should have DIFFERENT UUIDs"
