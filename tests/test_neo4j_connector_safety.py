"""Credential material must never reach a log, an error, or an artifact.

Both repositories are public. `_safe_neo4j_url()` redacts our own message, but
every caller used to append py2neo's raw exception text, whose ConnectionProfile
repr re-emits the password fragment the redactor exists to avoid. See issue #62.
"""

import os
import sys

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


# --- Boundary redaction -----------------------------------------------------
#
# An adversarial review showed that per-call-site gating is not enough: an
# ungated exc_info=True in reaction_generator, a chained __cause__, and two
# entrypoints with no handler at all each still emitted
# ConnectionProfile('bolt://ssSECRET@host') into debug_log.txt or stderr.
# Redaction now happens where text is emitted, so these test that boundary.

import logging

from src.credential_redaction import CredentialRedactingFilter, scrub


@pytest.mark.parametrize(
    "text",
    [
        "Cannot open connection to ConnectionProfile('bolt://ssSECRET@localhost:7687')",
        "Cannot connect to IPv4Address(('ssSECRET@localhost', 7687))",
        "connecting to bolt://neo4j:SECRET@localhost:7687",
        "profile ConnectionProfile(\"bolt://neo4j:SECRET@h\") failed",
    ],
)
def test_scrub_removes_credential_material(text):
    assert "SECRET" not in scrub(text)


def test_scrub_leaves_ordinary_text_alone():
    msg = "Reaction R-HSA-69620 has empty outputs, skipping"
    assert scrub(msg) == msg


def test_scrub_replaces_the_password_without_discarding_the_message(monkeypatch):
    """A common-word password must not blank out the whole line."""
    monkeypatch.setenv("NEO4J_PASSWORD", "reactome")
    out = scrub("Cannot find /opt/reactome/data/graph.db")
    assert "reactome" not in out
    assert "Cannot find" in out and "graph.db" in out


def test_scrub_ignores_trivially_short_passwords(monkeypatch):
    """Substituting a 1-3 char value would corrupt every line for no benefit."""
    monkeypatch.setenv("NEO4J_PASSWORD", "ab")
    assert scrub("a table of abbreviations") == "a table of abbreviations"


def test_logging_filter_scrubs_message_and_traceback(caplog):
    """The path that defeated per-call-site gating: an ungated exc_info."""
    record_filter = CredentialRedactingFilter()
    try:
        raise RuntimeError(
            "Cannot open connection to ConnectionProfile('bolt://ssSECRET@h')"
        )
    except RuntimeError:
        record = logging.LogRecord(
            "t", logging.ERROR, __file__, 1,
            "failed talking to bolt://neo4j:ssSECRET@h", (), sys.exc_info(),
        )
    assert record_filter.filter(record) is True
    assert "SECRET" not in record.getMessage()
    assert "SECRET" not in (record.exc_text or "")
    # exc_info must be cleared, or the handler re-renders the original.
    assert record.exc_info is None


def test_safe_exception_scrubs_rather_than_withholds_a_password_substring(monkeypatch):
    monkeypatch.setenv("NEO4J_PASSWORD", "reactome")
    rendered = nc._safe_exception(RuntimeError("KeyError: 'reactome_release'"))
    assert "reactome" not in rendered
    assert "withheld" not in rendered
