"""Tests for utility functions that were previously untested."""

import pytest
import pandas as pd
import numpy as np
from typing import Any
import sys
from pathlib import Path
from unittest.mock import patch

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import functions to test
from src.reaction_generator import is_valid_uuid
from src.logic_network_generator import (
    _get_reactome_id_from_hash,
    _get_hash_for_reaction,
    _get_non_null_values
)


class TestIsValidUUID:
    """Test the is_valid_uuid function."""

    def test_valid_64_char_string(self):
        """Valid UUID is 64-character string."""
        valid_uuid = "a" * 64
        assert is_valid_uuid(valid_uuid) is True

    def test_invalid_short_string(self):
        """String shorter than 64 characters is invalid."""
        short_uuid = "a" * 63
        assert is_valid_uuid(short_uuid) is False

    def test_invalid_long_string(self):
        """String longer than 64 characters is invalid."""
        long_uuid = "a" * 65
        assert is_valid_uuid(long_uuid) is False

    def test_empty_string(self):
        """Empty string is invalid."""
        assert is_valid_uuid("") is False

    def test_none_value(self):
        """None value should return False, not crash."""
        assert is_valid_uuid(None) is False

    def test_integer_value(self):
        """Integer value should return False, not crash."""
        assert is_valid_uuid(12345) is False

    def test_list_value(self):
        """List value should return False, not crash."""
        assert is_valid_uuid([]) is False

    def test_dict_value(self):
        """Dict value should return False, not crash."""
        assert is_valid_uuid({}) is False

    def test_actual_hash_format(self):
        """Test with actual SHA256-like hash."""
        sha256_hash = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        assert is_valid_uuid(sha256_hash) is True

    def test_hex_string_wrong_length(self):
        """Hex string with wrong length is invalid."""
        hex_string = "abc123"
        assert is_valid_uuid(hex_string) is False


class TestGetReactomeIdFromHash:
    """Test _get_reactome_id_from_hash function."""

    def test_successful_lookup(self):
        """Test successful hash lookup."""
        df = pd.DataFrame({
            "uid": ["hash1", "hash2", "hash3"],
            "reactome_id": ["R-HSA-100", "R-HSA-200", "R-HSA-300"]
        })
        result = _get_reactome_id_from_hash(df, "hash2")
        assert result == "R-HSA-200"

    def test_first_hash_lookup(self):
        """Test lookup of first hash."""
        df = pd.DataFrame({
            "uid": ["hash1", "hash2"],
            "reactome_id": ["R-HSA-100", "R-HSA-200"]
        })
        result = _get_reactome_id_from_hash(df, "hash1")
        assert result == "R-HSA-100"

    def test_last_hash_lookup(self):
        """Test lookup of last hash."""
        df = pd.DataFrame({
            "uid": ["hash1", "hash2", "hash3"],
            "reactome_id": ["R-HSA-100", "R-HSA-200", "R-HSA-300"]
        })
        result = _get_reactome_id_from_hash(df, "hash3")
        assert result == "R-HSA-300"

    def test_missing_hash_raises_error(self):
        """Missing hash should raise IndexError."""
        df = pd.DataFrame({
            "uid": ["hash1", "hash2"],
            "reactome_id": ["R-HSA-100", "R-HSA-200"]
        })
        with pytest.raises(IndexError):
            _get_reactome_id_from_hash(df, "nonexistent")

    def test_empty_dataframe_raises_error(self):
        """Empty DataFrame should raise IndexError."""
        df = pd.DataFrame({
            "uid": [],
            "reactome_id": []
        })
        with pytest.raises(IndexError):
            _get_reactome_id_from_hash(df, "any_hash")

    def test_duplicate_hashes_returns_first(self):
        """When duplicate hashes exist, returns first match."""
        df = pd.DataFrame({
            "uid": ["hash1", "hash1", "hash2"],
            "reactome_id": ["R-HSA-100", "R-HSA-999", "R-HSA-200"]
        })
        result = _get_reactome_id_from_hash(df, "hash1")
        # Should return first match
        assert result == "R-HSA-100"


class TestGetHashForReaction:
    """Test _get_hash_for_reaction function."""

    def test_successful_input_hash_lookup(self):
        """Test successful lookup of input hash."""
        df = pd.DataFrame({
            "uid": ["uid1", "uid2"],
            "input_hash": ["hash_in1", "hash_in2"],
            "output_hash": ["hash_out1", "hash_out2"]
        })
        result = _get_hash_for_reaction(df, "uid2", "input_hash")
        assert result == "hash_in2"

    def test_successful_output_hash_lookup(self):
        """Test successful lookup of output hash."""
        df = pd.DataFrame({
            "uid": ["uid1", "uid2"],
            "input_hash": ["hash_in1", "hash_in2"],
            "output_hash": ["hash_out1", "hash_out2"]
        })
        result = _get_hash_for_reaction(df, "uid1", "output_hash")
        assert result == "hash_out1"

    def test_missing_uid_raises_error(self):
        """Missing UID should raise IndexError."""
        df = pd.DataFrame({
            "uid": ["uid1", "uid2"],
            "input_hash": ["hash1", "hash2"]
        })
        with pytest.raises(IndexError):
            _get_hash_for_reaction(df, "nonexistent", "input_hash")

    def test_empty_dataframe_raises_error(self):
        """Empty DataFrame should raise IndexError."""
        df = pd.DataFrame({
            "uid": [],
            "input_hash": []
        })
        with pytest.raises(IndexError):
            _get_hash_for_reaction(df, "any_uid", "input_hash")


class TestGetNonNullValues:
    """Test _get_non_null_values function."""

    def test_all_non_null_values(self):
        """All non-null values are returned."""
        df = pd.DataFrame({"col": [1, 2, 3]})
        result = _get_non_null_values(df, "col")
        assert result == [1, 2, 3]

    def test_removes_none_values(self):
        """None values are filtered out."""
        df = pd.DataFrame({"col": [1, None, 2, None, 3]})
        result = _get_non_null_values(df, "col")
        assert result == [1, 2, 3]

    def test_removes_nan_values(self):
        """NaN values are filtered out."""
        df = pd.DataFrame({"col": [1, np.nan, 2, np.nan, 3]})
        result = _get_non_null_values(df, "col")
        assert result == [1, 2, 3]

    def test_empty_dataframe(self):
        """Empty DataFrame returns empty list."""
        df = pd.DataFrame({"col": []})
        result = _get_non_null_values(df, "col")
        assert result == []

    def test_all_null_values(self):
        """Column of all null values returns empty list."""
        df = pd.DataFrame({"col": [None, np.nan, None]})
        result = _get_non_null_values(df, "col")
        assert result == []

    def test_preserves_order(self):
        """Non-null values maintain their original order."""
        df = pd.DataFrame({"col": [3, None, 1, None, 2]})
        result = _get_non_null_values(df, "col")
        assert result == [3, 1, 2]

    def test_handles_zero(self):
        """Zero is not treated as null."""
        df = pd.DataFrame({"col": [0, None, 1, None, 2]})
        result = _get_non_null_values(df, "col")
        assert result == [0, 1, 2]

    def test_handles_empty_string(self):
        """Empty string is not treated as null."""
        df = pd.DataFrame({"col": ["", None, "a", None, "b"]})
        result = _get_non_null_values(df, "col")
        assert result == ["", "a", "b"]

    def test_handles_false(self):
        """False is not treated as null."""
        df = pd.DataFrame({"col": [False, None, True, None, False]})
        result = _get_non_null_values(df, "col")
        assert result == [False, True, False]


class TestDataFrameEdgeCases:
    """Test edge cases with DataFrames."""

    def test_dataframe_with_missing_columns(self):
        """DataFrame missing expected columns should raise KeyError."""
        df = pd.DataFrame({
            "wrong_column": ["value1", "value2"]
        })
        with pytest.raises(KeyError):
            _get_reactome_id_from_hash(df, "hash1")

    def test_dataframe_with_null_values_in_uid(self):
        """DataFrame with null UIDs should not match."""
        import numpy as np
        df = pd.DataFrame({
            "uid": ["hash1", np.nan, "hash3"],
            "reactome_id": ["R-HSA-100", "R-HSA-200", "R-HSA-300"]
        })
        with pytest.raises(IndexError):
            # np.nan != np.nan, so this should not match
            _get_reactome_id_from_hash(df, np.nan)

    def test_dataframe_with_duplicate_columns(self):
        """DataFrame can have duplicate column names (pandas allows this)."""
        # This is more of a pandas quirk test
        df = pd.DataFrame({
            "uid": ["hash1", "hash2"],
            "reactome_id": ["R-HSA-100", "R-HSA-200"]
        })
        # Just verify it works normally
        result = _get_reactome_id_from_hash(df, "hash1")
        assert result == "R-HSA-100"


class TestTypeConversions:
    """Test type conversion edge cases."""

    def test_stable_id_returned_as_string(self):
        """Reactome stable ID should be returned as string."""
        df = pd.DataFrame({
            "uid": ["hash1"],
            "reactome_id": ["R-HSA-100"]
        })
        result = _get_reactome_id_from_hash(df, "hash1")
        assert isinstance(result, str)
        assert result == "R-HSA-100"

    def test_string_uid_comparison(self):
        """UID comparison should work with strings."""
        df = pd.DataFrame({
            "uid": ["hash1", "hash2"],
            "reactome_id": ["R-HSA-100", "R-HSA-200"]
        })
        result = _get_reactome_id_from_hash(df, "hash1")
        assert result == "R-HSA-100"

    def test_numeric_string_uid(self):
        """Numeric string UIDs should work."""
        df = pd.DataFrame({
            "uid": ["123", "456"],
            "reactome_id": ["R-HSA-100", "R-HSA-200"]
        })
        result = _get_reactome_id_from_hash(df, "456")
        assert result == "R-HSA-200"
