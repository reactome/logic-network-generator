"""Credential material must never reach a log, an error, or an artifact.

Both repositories are public. `_safe_neo4j_url()` redacts our own message, but
every caller used to append py2neo's raw exception text, whose ConnectionProfile
repr re-emits the password fragment the redactor exists to avoid. See issue #62.
"""

import os

import pytest

import src.neo4j_connector as nc


@pytest.mark.parametrize(
    "url,expected",
    [
        ("bolt://neo4j:hunter2@localhost:7687", "bolt://localhost:7687"),
        # An unencoded "@" in the password: ConnectionProfile.uri emitted the
        # tail of the password as if it were the host.
        ("bolt://neo4j:p@ssSECRET@localhost:7687", "bolt://localhost:7687"),
        # An unencoded "/" makes the authority boundary undecidable, so we
        # disclose nothing rather than risk a fragment.
        ("bolt://neo4j:pw/slashSECRET@localhost:7687", "bolt://<redacted>"),
        ("bolt://neo4j:hunter2@localhost:7687#TOPSECRET", "bolt://localhost:7687"),
        ("bolt://localhost:7687?password=hunter2", "bolt://localhost:7687"),
        ("bolt://[::1]:7687", "bolt://[::1]:7687"),
        ("not-a-url", "<neo4j-url>"),
    ],
)
def test_safe_url_never_emits_credential_material(url, expected, monkeypatch):
    monkeypatch.setenv("NEO4J_URL", url)
    got = nc._safe_neo4j_url()
    assert got == expected
    assert "SECRET" not in got
    assert "hunter2" not in got


@pytest.mark.parametrize(
    "text",
    [
        "Cannot open connection to ConnectionProfile('bolt://ssSECRET@host')",
        "Cannot connect to IPv4Address(('ssSECRET@localhost', 9999))",
        "auth failure for bolt://neo4j:SECRET@host",
    ],
)
def test_credential_bearing_exception_text_is_withheld(text):
    rendered = nc._safe_exception(Exception(text))
    assert "SECRET" not in rendered
    assert "withheld" in rendered


def test_exception_text_matching_the_password_is_withheld(monkeypatch):
    monkeypatch.setenv("NEO4J_PASSWORD", "hunter2")
    assert "hunter2" not in nc._safe_exception(Exception("auth failed for hunter2"))


def test_benign_exception_keeps_its_message():
    """Over-redaction would make every failure undiagnosable."""
    rendered = nc._safe_exception(ValueError("column 'foo' not found in results"))
    assert rendered == "ValueError: column 'foo' not found in results"


def test_tracebacks_are_off_by_default(monkeypatch):
    """A py2neo traceback carries the same ConnectionProfile the message did."""
    monkeypatch.delenv("LNG_DEBUG_TRACEBACKS", raising=False)
    assert nc._traceback_kwargs() == {"exc_info": False}
    monkeypatch.setenv("LNG_DEBUG_TRACEBACKS", "1")
    assert nc._traceback_kwargs() == {"exc_info": True}
