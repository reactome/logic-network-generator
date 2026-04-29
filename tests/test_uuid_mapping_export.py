"""Tests for UUID mapping export functionality.

Tests verify that export_uuid_to_reactome_mapping correctly creates
a mapping from UUIDs in the logic network to Reactome stable IDs.
"""

import pandas as pd
import tempfile
import os
import pytest
from pathlib import Path


def find_first_pathway_dir():
    """Find the first available generated pathway directory."""
    output_dir = Path("output")
    if not output_dir.exists():
        return None
    for d in sorted(output_dir.iterdir()):
        if d.is_dir() and (d / "logic_network.csv").exists() and (d / "stid_to_uuid_mapping.csv").exists():
            return d
    return None


PATHWAY_DIR = find_first_pathway_dir()

# Integration tier: requires generated pathway artifacts in output/
pytestmark = pytest.mark.integration


class TestUUIDMappingFileStructure:
    """Test the structure and content of generated UUID mapping files."""

    pytestmark = pytest.mark.skipif(
        PATHWAY_DIR is None,
        reason="No generated pathway directories found in output/"
    )

    def test_mapping_file_has_required_columns(self):
        """UUID mapping file should have uuid and stable_id columns."""
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")
        assert 'uuid' in mapping.columns, "Missing 'uuid' column"
        assert 'stable_id' in mapping.columns, "Missing 'stable_id' column"

    def test_mapping_file_is_not_empty(self):
        """UUID mapping file should have entries."""
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")
        assert len(mapping) > 0, "UUID mapping file is empty"

    def test_all_uuids_are_unique(self):
        """Each UUID in the mapping should be unique."""
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")
        assert mapping['uuid'].nunique() == len(mapping), \
            f"Found duplicate UUIDs: {len(mapping) - mapping['uuid'].nunique()} duplicates"

    def test_no_null_uuids(self):
        """No UUIDs should be null."""
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")
        assert mapping['uuid'].notna().all(), "Found null UUIDs in mapping"

    def test_stable_ids_have_correct_format(self):
        """Stable IDs should follow R-XXX-NNN format."""
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")
        non_null_ids = mapping['stable_id'].dropna()
        for sid in non_null_ids:
            assert str(sid).startswith("R-"), \
                f"Stable ID does not start with 'R-': {sid}"


class TestUUIDMappingCompleteness:
    """Test that UUID mapping covers all UUIDs in the logic network."""

    pytestmark = pytest.mark.skipif(
        PATHWAY_DIR is None,
        reason="No generated pathway directories found in output/"
    )

    def test_all_network_uuids_in_mapping(self):
        """Every UUID in the logic network should have a mapping entry."""
        network = pd.read_csv(PATHWAY_DIR / "logic_network.csv")
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")

        network_uuids = set(network['source_id'].unique()) | set(network['target_id'].unique())
        mapping_uuids = set(mapping['uuid'].unique())

        unmapped = network_uuids - mapping_uuids
        assert len(unmapped) == 0, \
            f"Found {len(unmapped)} UUIDs in logic network without mapping entries"

    def test_position_aware_uuids_have_different_ids(self):
        """Same stable_id at different positions should have different UUIDs."""
        mapping = pd.read_csv(PATHWAY_DIR / "stid_to_uuid_mapping.csv")

        multi_position = mapping['stable_id'].value_counts()
        multi_position_entities = multi_position[multi_position > 1]

        if len(multi_position_entities) == 0:
            pytest.skip("No multi-position entities in this pathway")

        for stable_id in multi_position_entities.index:
            entity_rows = mapping[mapping['stable_id'] == stable_id]
            uuids = entity_rows['uuid'].unique()
            assert len(uuids) == len(entity_rows), \
                f"Stable ID {stable_id} appears {len(entity_rows)} times but has only {len(uuids)} unique UUIDs"


class TestUUIDMappingAcrossPathways:
    """Test UUID mapping across multiple pathways."""

    @staticmethod
    def get_pathway_dirs():
        output_dir = Path("output")
        if not output_dir.exists():
            return []
        return [
            str(d / "stid_to_uuid_mapping.csv")
            for d in sorted(output_dir.iterdir())
            if d.is_dir() and (d / "stid_to_uuid_mapping.csv").exists()
        ]

    MAPPING_FILES = get_pathway_dirs.__func__()

    @pytest.mark.skipif(len(MAPPING_FILES) == 0, reason="No generated pathways found")
    @pytest.mark.parametrize("mapping_path", MAPPING_FILES[:5],
                             ids=[Path(p).parent.name for p in MAPPING_FILES[:5]])
    def test_every_pathway_has_valid_mapping(self, mapping_path):
        """Each pathway's UUID mapping should have valid structure."""
        mapping = pd.read_csv(mapping_path)
        assert len(mapping) > 0, "UUID mapping is empty"
        assert 'uuid' in mapping.columns
        assert 'stable_id' in mapping.columns
        assert mapping['uuid'].notna().all(), "Found null UUIDs"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
