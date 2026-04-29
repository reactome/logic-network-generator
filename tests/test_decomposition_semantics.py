"""Lockdown tests for the decomposition contract in src/reaction_generator.py.

These tests codify the rules a future reader (human or LLM) needs to obey
before "fixing" decomposition. The rules are biologically grounded — see
docs/DESIGN_DECISIONS.md, "Complex vs EntitySet" — and a previous round of
this work nearly broke them by collapsing the two into one concept. If you
are about to change one of these tests, read that doc first.

The contract:

- EntitySet → flat set of alternatives. NEVER a cartesian product.
- Simple Complex (no EntitySet inside) → returned intact, not decomposed.
- Complex containing EntitySet → cartesian-product UIDs over combinations.
- Ubiquitin EntitySets → returned intact (combinatorial-explosion guard).
- Simple entity → returned intact.

These rules govern `decomposed_uid_mapping` only — the plumbing for
cross-reaction matching. The final `logic_network.csv` has its own,
different rules around root/terminal complex decomposition.
"""

from unittest.mock import patch

import pytest

import src.reaction_generator as rg


@pytest.fixture(autouse=True)
def reset_module_state():
    """Each test runs against a fresh decomposition store and caches."""
    rg._store.clear()
    rg.reference_entity_dict.clear()
    rg._complex_contains_set_cache.clear()
    rg._direct_component_stoichiometry.clear()
    yield


def _store_df():
    """Materialize the decomposition store as a DataFrame for assertions.

    The store now holds rows in a list with dict indexes; tests pull a
    DataFrame view at assertion time so the existing filter expressions
    still read naturally.
    """
    return rg._store.to_dataframe()


def _label_map(mapping):
    """Build a get_labels stub from {entity_id: [labels]}."""
    def _f(entity_id):
        return mapping.get(entity_id, ["GenomeEncodedEntity"])
    return _f


def _components_map(mapping):
    """Build a get_complex_components stub from {complex_id: {member: stoich}}."""
    def _f(entity_id):
        return mapping.get(entity_id, {})
    return _f


def _members_map(mapping):
    """Build a get_set_members stub from {set_id: [members]}."""
    def _f(entity_id):
        return mapping.get(entity_id, [])
    return _f


def _ref_entity_map(mapping):
    def _f(entity_id):
        return mapping.get(entity_id)
    return _f


class TestEntitySetIsFlatAlternatives:
    """EntitySet members are alternatives — flat set, never cartesian."""

    def test_entityset_returns_flat_set_of_simple_members(self):
        labels = _label_map({"S": ["EntitySet"]})
        members = _members_map({"S": ["A", "B", "C"]})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_set_members", members), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            result = rg.break_apart_entity("S")

        assert result == {"A", "B", "C"}, (
            "EntitySet must return alternatives flat. "
            "If this fails returning a cartesian product or UIDs, the "
            "matcher will produce false links between alternatives."
        )

    def test_entityset_writes_provenance_rows_with_source_entity_id(self):
        labels = _label_map({"S": ["EntitySet"]})
        members = _members_map({"S": ["A", "B"]})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_set_members", members), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            rg.break_apart_entity("S")

        df = _store_df()
        prov = df[df["source_entity_id"] == "S"]
        assert set(prov["component_id"]) == {"A", "B"}, (
            "EntitySet decomposition must record provenance for each leaf so "
            "downstream queries can answer 'which set does leaf X belong to?'"
        )


class TestSimpleComplexIsAtomic:
    """Simple complex (no EntitySet) is one species; matcher must NOT see through."""

    def test_simple_complex_returns_itself(self):
        labels = _label_map({"C": ["Complex"], "A": ["GenomeEncodedEntity"], "B": ["GenomeEncodedEntity"]})
        components = _components_map({"C": {"A": 1, "B": 1}})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_complex_components", components), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            result = rg.break_apart_entity("C")

        assert result == {"C"}, (
            "Simple complex MUST return itself intact. Decomposing it lets the "
            "matcher link unrelated species (e.g. complex A:B → free A) "
            "through component overlap, which is biologically wrong."
        )

    def test_simple_complex_writes_no_decomposition_rows(self):
        labels = _label_map({"C": ["Complex"], "A": ["GenomeEncodedEntity"], "B": ["GenomeEncodedEntity"]})
        components = _components_map({"C": {"A": 1, "B": 1}})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_complex_components", components), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            rg.break_apart_entity("C")

        df = _store_df()
        rows_for_c = df[df["reactome_id"] == "C"]
        assert len(rows_for_c) == 0, (
            "Atomic complex must not write decomposition rows."
        )


class TestComplexWithEntitySetIsCartesian:
    """A Complex that contains an EntitySet IS decomposed — alternatives need expansion."""

    def test_complex_with_entityset_returns_uids(self):
        labels = _label_map({
            "C": ["Complex"],
            "S": ["EntitySet"],
            "P": ["GenomeEncodedEntity"],
            "A": ["GenomeEncodedEntity"],
            "B": ["GenomeEncodedEntity"],
        })
        components = _components_map({"C": {"S": 1, "P": 1}})
        members = _members_map({"S": ["A", "B"]})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_complex_components", components), \
             patch.object(rg, "get_set_members", members), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            result = rg.break_apart_entity("C")

        # Cartesian product over {A, B} × {P} = two combinations: {A, P} and {B, P}.
        # Each combination gets a UID (sha256 hash, 64 chars).
        assert all(len(uid) == 64 for uid in result), (
            "Complex-with-EntitySet must return UIDs for combinations, not raw IDs"
        )
        assert len(result) == 2, (
            f"Cartesian product of {{A,B}} × {{P}} should yield 2 combinations, got {len(result)}"
        )

    def test_complex_with_entityset_writes_rows_per_combination(self):
        labels = _label_map({
            "C": ["Complex"],
            "S": ["EntitySet"],
            "P": ["GenomeEncodedEntity"],
            "A": ["GenomeEncodedEntity"],
            "B": ["GenomeEncodedEntity"],
        })
        components = _components_map({"C": {"S": 1, "P": 1}})
        members = _members_map({"S": ["A", "B"]})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_complex_components", components), \
             patch.object(rg, "get_set_members", members), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            rg.break_apart_entity("C")

        # Each combination has 2 components → 4 rows total under reactome_id=C
        df = _store_df()
        rows_for_c = df[df["reactome_id"] == "C"]
        assert len(rows_for_c) == 4
        assert set(rows_for_c["component_id"]) == {"A", "B", "P"}


class TestUbiquitinIsAtomic:
    """Ubiquitin EntitySets are skipped to avoid combinatorial explosion."""

    def test_ubiquitin_entityset_returned_intact(self):
        # R-HSA-113595 is one of the registered ubiquitin EntitySet IDs.
        labels = _label_map({"R-HSA-113595": ["EntitySet"]})
        with patch.object(rg, "get_labels", labels):
            result = rg.break_apart_entity("R-HSA-113595")
        assert result == {"R-HSA-113595"}, (
            "Ubiquitin EntitySets must NOT be decomposed; their members are "
            "near-identical 76-aa proteins and decomposing them causes a "
            "combinatorial explosion with no biological insight gained."
        )


class TestSimpleEntityIsAtomic:
    """Simple entities (proteins, small molecules) are returned as-is."""

    @pytest.mark.parametrize("label", [
        "GenomeEncodedEntity",
        "EntityWithAccessionedSequence",
        "SimpleEntity",
        "ChemicalDrug",
        "Polymer",
        "OtherEntity",
        "Cell",
        "Drug",
    ])
    def test_simple_entity_returns_itself(self, label):
        labels = _label_map({"X": [label]})
        with patch.object(rg, "get_labels", labels):
            result = rg.break_apart_entity("X")
        assert result == {"X"}


class TestCrossCallStability:
    """Same entity decomposed twice returns the same leaves (cache short-circuit)."""

    def test_repeated_entityset_decomposition_is_stable(self):
        labels = _label_map({"S": ["EntitySet"]})
        members = _members_map({"S": ["A", "B"]})
        with patch.object(rg, "get_labels", labels), \
             patch.object(rg, "get_set_members", members), \
             patch.object(rg, "get_reference_entity_id", _ref_entity_map({})):
            first = rg.break_apart_entity("S")
            second = rg.break_apart_entity("S")
        assert first == second
        # Provenance rows should NOT be duplicated by the second call.
        df = _store_df()
        prov = df[df["source_entity_id"] == "S"]
        assert len(prov) == 2, (
            "Cache short-circuit must prevent duplicate provenance rows on re-entry."
        )
