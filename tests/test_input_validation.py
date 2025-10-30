"""Tests for input validation in create_pathway_logic_network."""

import pytest
import pandas as pd
import sys
from unittest.mock import patch

sys.path.insert(0, '/home/awright/gitroot/logic-network-generator')

# Mock py2neo.Graph to avoid Neo4j connection during import
with patch('py2neo.Graph'):
    from src.logic_network_generator import create_pathway_logic_network


class TestInputValidation:
    """Test that create_pathway_logic_network validates its inputs properly."""

    def test_rejects_empty_decomposed_uid_mapping(self):
        """Should raise ValueError if decomposed_uid_mapping is empty."""
        empty_mapping = pd.DataFrame()
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="decomposed_uid_mapping cannot be empty"):
            create_pathway_logic_network(empty_mapping, valid_connections, valid_matches)

    def test_rejects_decomposed_uid_mapping_missing_uid_column(self):
        """Should raise ValueError if decomposed_uid_mapping is missing 'uid' column."""
        invalid_mapping = pd.DataFrame({
            # Missing 'uid' column
            'reactome_id': [1, 2],
            'input_or_output_reactome_id': [10, 20]
        })
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="missing required columns.*uid"):
            create_pathway_logic_network(invalid_mapping, valid_connections, valid_matches)

    def test_rejects_decomposed_uid_mapping_missing_reactome_id_column(self):
        """Should raise ValueError if decomposed_uid_mapping is missing 'reactome_id' column."""
        invalid_mapping = pd.DataFrame({
            'uid': ['hash1', 'hash2'],
            # Missing 'reactome_id' column
            'input_or_output_reactome_id': [10, 20]
        })
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="missing required columns.*reactome_id"):
            create_pathway_logic_network(invalid_mapping, valid_connections, valid_matches)

    def test_rejects_decomposed_uid_mapping_missing_input_or_output_reactome_id_column(self):
        """Should raise ValueError if missing 'input_or_output_reactome_id' column."""
        invalid_mapping = pd.DataFrame({
            'uid': ['hash1', 'hash2'],
            'reactome_id': [1, 2],
            # Missing 'input_or_output_reactome_id' column
        })
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="missing required columns.*input_or_output_reactome_id"):
            create_pathway_logic_network(invalid_mapping, valid_connections, valid_matches)

    def test_rejects_empty_reaction_connections(self):
        """Should raise ValueError if reaction_connections is empty."""
        valid_mapping = pd.DataFrame({
            'uid': ['hash1', 'hash2'],
            'reactome_id': [1, 2],
            'input_or_output_reactome_id': [10, 20],
            'component_id': [0, 0],
            'component_id_or_reference_entity_id': [0, 0],
            'input_or_output_uid': [None, None]
        })
        empty_connections = pd.DataFrame()
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="reaction_connections cannot be empty"):
            create_pathway_logic_network(valid_mapping, empty_connections, valid_matches)

    def test_rejects_reaction_connections_missing_preceding_reaction_id(self):
        """Should raise ValueError if reaction_connections is missing 'preceding_reaction_id'."""
        valid_mapping = pd.DataFrame({
            'uid': ['hash1', 'hash2'],
            'reactome_id': [1, 2],
            'input_or_output_reactome_id': [10, 20],
            'component_id': [0, 0],
            'component_id_or_reference_entity_id': [0, 0],
            'input_or_output_uid': [None, None]
        })
        invalid_connections = pd.DataFrame({
            # Missing 'preceding_reaction_id'
            'following_reaction_id': [2, 3]
        })
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="missing required columns.*preceding_reaction_id"):
            create_pathway_logic_network(valid_mapping, invalid_connections, valid_matches)

    def test_rejects_empty_best_matches(self):
        """Should raise ValueError if best_matches is empty DataFrame."""
        valid_mapping = pd.DataFrame({
            'uid': ['hash1', 'hash2'],
            'reactome_id': [1, 2],
            'input_or_output_reactome_id': [10, 20],
            'component_id': [0, 0],
            'component_id_or_reference_entity_id': [0, 0],
            'input_or_output_uid': [None, None]
        })
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        empty_matches = pd.DataFrame()

        with pytest.raises(ValueError, match="best_matches cannot be empty"):
            create_pathway_logic_network(valid_mapping, valid_connections, empty_matches)

    def test_rejects_best_matches_missing_incomming_column(self):
        """Should raise ValueError if best_matches is missing 'incomming' column."""
        valid_mapping = pd.DataFrame({
            'uid': ['hash1', 'hash2'],
            'reactome_id': [1, 2],
            'input_or_output_reactome_id': [10, 20],
            'component_id': [0, 0],
            'component_id_or_reference_entity_id': [0, 0],
            'input_or_output_uid': [None, None]
        })
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        invalid_matches = pd.DataFrame({
            # Missing 'incomming' column
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError, match="missing required columns.*incomming"):
            create_pathway_logic_network(valid_mapping, valid_connections, invalid_matches)

    def test_error_message_shows_available_columns(self):
        """Error messages should show what columns are actually available."""
        invalid_mapping = pd.DataFrame({
            'wrong_column': [1, 2],
            'another_wrong_column': [3, 4]
        })
        valid_connections = pd.DataFrame({
            'preceding_reaction_id': [1, 2],
            'following_reaction_id': [2, 3]
        })
        valid_matches = pd.DataFrame({
            'incomming': ['hash1', 'hash2'],
            'outgoing': ['hash3', 'hash4']
        })

        with pytest.raises(ValueError) as exc_info:
            create_pathway_logic_network(invalid_mapping, valid_connections, valid_matches)

        error_msg = str(exc_info.value)
        assert "Available columns:" in error_msg
        assert "wrong_column" in error_msg
        assert "another_wrong_column" in error_msg
