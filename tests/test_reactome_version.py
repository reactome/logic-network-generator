"""Reactome version sentinel.

This isn't a correctness test — it's a recorder. It runs a single Cypher
query that asks Neo4j what release it has loaded, then asserts a sane
range and prints the version into the test output. Purpose:

  - Every database-tier CI run logs which Reactome the tests ran against.
    When a teammate sees a regression, they can tell whether it correlates
    with a Reactome release bump.
  - Pinning a docker-compose image is a guess; this confirms what's
    actually loaded.
  - When Reactome ships v97 (and we bump docker-compose), this test still
    passes — it's not version-specific. It just records the new number.

If this test fails because the version is unexpectedly old or absent,
docker-compose was probably never started or is pointing at the wrong
image.
"""

import pytest


@pytest.mark.database
def test_reactome_version_loaded(capsys):
    """Read the loaded Reactome release and print it into the test log."""
    from src.neo4j_connector import get_graph

    graph = get_graph()
    rows = graph.run(
        """
        MATCH (info:DBInfo)
        RETURN info.name AS name,
               toInteger(coalesce(info.releaseNumber, info.version)) AS version
        LIMIT 1
        """
    ).data()

    assert rows, (
        "Neo4j has no DBInfo node. Either the database isn't loaded with "
        "Reactome data, or it's a non-Reactome graph."
    )

    version = rows[0]["version"]
    name = rows[0]["name"]
    assert isinstance(version, int) and version >= 80, (
        f"Reactome version {version!r} is implausibly old; check that "
        f"docker-compose is pointing at a recent Release tag."
    )

    # Print into pytest's capture buffer so CI logs show what we tested.
    with capsys.disabled():
        print(f"\n[reactome-version] tests ran against {name} v{version}")
