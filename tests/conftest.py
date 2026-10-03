"""Shared fixtures.

Removed LNG_* flags are rejected at runtime (see `_REMOVED_ENV`), so a stale
value exported in a developer's shell or a CI runner would turn every test that
builds a network into a spurious failure. Tests should exercise the code, not
the ambient environment, so clear those names for the whole suite. The guard
itself is tested explicitly in tests/test_removed_env_flags.py.
"""

import pytest

from src.logic_network_generator import _REMOVED_ENV


@pytest.fixture(autouse=True)
def _clear_removed_env_flags(monkeypatch):
    for name in _REMOVED_ENV:
        monkeypatch.delenv(name, raising=False)
    yield


@pytest.fixture(autouse=True)
def _offline_neo4j_lookups(request, monkeypatch):
    """Unit tests run without Neo4j (CI has none). Two lookups used to fall back
    silently when the database was unreachable, and the unit tests leaned on
    that: offline they got the fallback, and on a developer machine they
    quietly queried the live database. The fallbacks now raise (code review
    2026-10-02), so the offline answer is stated here, explicitly. Tests marked
    `database` keep the real lookups."""
    if request.node.get_closest_marker("database"):
        yield
        return
    import src.neo4j_connector as nc
    import src.reaction_generator as rg
    monkeypatch.setattr(nc, "get_modifier_isoform_entity_set_ids",
                        lambda: set(rg._UBIQUITIN_ENTITY_SET_IDS))
    monkeypatch.setattr(nc, "get_pathway_participating_entities", lambda pathway_id: set())
    monkeypatch.setattr(rg, "_modifier_set_cache", None)
    yield
