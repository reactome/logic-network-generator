"""Redact credential material at the output boundary.

Both Reactome repositories are public, and the realistic disclosure path is
someone pasting a log line or a stack trace into an issue. Neo4j credentials
can be embedded in ``NEO4J_URL``, and py2neo re-emits them in its own
exception text via ``ConnectionProfile`` and ``IPv4Address`` reprs.

Redacting at each call site does not work, and an adversarial review proved
it: gating the fourteen ``exc_info=True`` calls inside ``neo4j_connector``
left the same secret reaching ``debug_log.txt`` through an ungated
``exc_info=True`` two frames up in ``reaction_generator``, and reaching stderr
entirely unhandled from ``scripts/validate_logic_network.py``. Chained
exceptions leak it again through ``__cause__`` even when ``str()`` is clean.

So redaction belongs where the text is emitted — a logging filter on every
handler, plus an excepthook for anything uncaught — which covers paths nobody
enumerated.

Scrubbing is structural rather than value-based wherever possible. Replacing
occurrences of the password itself is a poor primary defence here because the
documented dev passwords are ordinary words ("reactome", "test") that appear
in paths, keys and prose; value replacement is applied, but it substitutes the
token rather than discarding the whole message.
"""

from __future__ import annotations

import logging
import os
import re
import sys
import traceback
from typing import Any

REDACTED = "<redacted>"

# scheme://userinfo@host — the userinfo is credential material by definition.
_URL_USERINFO = re.compile(r"(?P<scheme>[A-Za-z][A-Za-z0-9+.\-]*://)[^\s/@]*@")

# py2neo's own reprs, which are how a password fragment escapes even when our
# own message is already redacted.
_PROFILE_REPR = re.compile(r"ConnectionProfile\((['\"]).*?\1\)")
_ADDRESS_REPR = re.compile(r"IPv4Address\(\(.*?\)\)")
_ADDRESS6_REPR = re.compile(r"IPv6Address\(\(.*?\)\)")


def scrub(text: str) -> str:
    """Remove credential material from arbitrary text.

    Structural first (URL userinfo, py2neo reprs), then the literal password
    as a backstop for paths that format it directly.
    """
    if not text:
        return text
    text = _URL_USERINFO.sub(lambda m: f"{m.group('scheme')}{REDACTED}@", text)
    text = _PROFILE_REPR.sub(f"ConnectionProfile({REDACTED})", text)
    text = _ADDRESS_REPR.sub(f"IPv4Address({REDACTED})", text)
    text = _ADDRESS6_REPR.sub(f"IPv6Address({REDACTED})", text)

    password = os.getenv("NEO4J_PASSWORD") or ""
    # Substituting a 1-3 character password would corrupt almost every line
    # for no security benefit; such a value is not a secret worth protecting.
    if len(password) >= 4 and password in text:
        text = text.replace(password, REDACTED)
    return text


class CredentialRedactingFilter(logging.Filter):
    """Scrub every log record, including its traceback, before it is emitted.

    Attached to handlers rather than loggers so it applies to records
    propagated from any module, including third-party ones.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        try:
            message = record.getMessage()
        except Exception:  # pragma: no cover - a broken record must still emit
            return True

        scrubbed = scrub(message)
        if scrubbed != message:
            record.msg = scrubbed
            record.args = ()

        if record.exc_info:
            # Render the traceback now and scrub it; leaving exc_info set
            # would let the handler re-render the unredacted original.
            rendered = "".join(traceback.format_exception(*record.exc_info))
            record.exc_text = scrub(rendered)
            record.exc_info = None
        elif record.exc_text:
            record.exc_text = scrub(record.exc_text)
        return True


def install_log_redaction(logger: logging.Logger | None = None) -> None:
    """Attach the filter to every handler on the root (or given) logger."""
    target = logger if logger is not None else logging.getLogger()
    for handler in target.handlers:
        if not any(isinstance(f, CredentialRedactingFilter) for f in handler.filters):
            handler.addFilter(CredentialRedactingFilter())


def install_excepthook() -> None:
    """Scrub uncaught tracebacks.

    `scripts/validate_logic_network.py` and `bin/validate-against-mpbiopath.py`
    construct a Graph with no handler at all, so a first-connection failure —
    the most likely failure — printed the full chain to stderr.
    """
    previous = sys.excepthook

    def _hook(exc_type: type, exc: BaseException, tb: Any) -> None:
        if previous is not sys.__excepthook__:
            previous(exc_type, exc, tb)
            return
        sys.stderr.write(scrub("".join(traceback.format_exception(exc_type, exc, tb))))

    sys.excepthook = _hook


def install() -> None:
    """Install both guards. Safe to call more than once."""
    install_log_redaction()
    install_excepthook()
