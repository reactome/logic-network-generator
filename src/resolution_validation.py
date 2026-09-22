"""Completeness checks for node_resolution.csv, in both directions.

Requirement: every node in the database must map to an LNG node, and back
again. "Perfectly" is not checkable, so it is expressed here as: no absence
goes undeclared.

Pure functions over already-loaded rows, deliberately free of the filesystem
and Neo4j, so the negative control in
``tests/test_resolution_negative_control.py`` can corrupt the inputs directly
and prove these checks FAIL. A completeness check that cannot fail is the
specific defect this project has shipped — ``_decomposed_ids`` once passed
11 of 11 while masking 18 dropped catalysts — so the negative control is the
load-bearing half of this module, not its polish.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, NamedTuple, Set


class Violation(NamedTuple):
    check: str
    subject: str
    detail: str


def check_reverse_completeness(network_node_uuids: Set[str],
                               resolution_rows: Iterable[dict]) -> List[Violation]:
    """Every node in the network must say what it stands for."""
    mapped = {str(r["uuid"]) for r in resolution_rows}
    return [
        Violation("reverse_completeness", uuid,
                  "node appears in the logic network but in no resolution row")
        for uuid in sorted(network_node_uuids - mapped)
    ]


def check_referential_integrity(network_node_uuids: Set[str],
                                resolution_rows: Iterable[dict]) -> List[Violation]:
    """Every resolution row must point at a node that exists.

    This is the shape of issue #67, where 100% of catalyst and regulator
    context rows named the undecomposed parent entity — a node absent from
    the network — and nothing noticed for as long as the export existed.
    """
    return [
        Violation("referential_integrity", str(r["uuid"]),
                  f"resolution row for {r['stable_id']} names a node absent "
                  f"from the logic network")
        for r in resolution_rows
        if str(r["uuid"]) not in network_node_uuids
    ]


def check_forward_completeness(participating_entities: Set[str],
                               resolution_rows: Iterable[dict],
                               exclusion_rows: Iterable[dict]) -> List[Violation]:
    """Every participating entity resolves, or is excluded with a reason."""
    resolved = {str(r["stable_id"]) for r in resolution_rows}
    excluded = {str(r["stable_id"]) for r in exclusion_rows}
    return [
        Violation("forward_completeness", stid,
                  "entity participates in a reaction but has no resolution "
                  "row and no exclusion")
        for stid in sorted(participating_entities - resolved - excluded)
    ]


def check_exclusions_have_reasons(exclusion_rows: Iterable[dict]) -> List[Violation]:
    """An exclusion without a reason is a silent absence wearing a hat.

    The reason is load-bearing rather than bureaucratic: the same list holds
    deliberately atomic modifier sets (a design decision) and entities that
    should have been generated (a bug), and only this column separates them.
    """
    violations = []
    for row in exclusion_rows:
        reason = str(row.get("reason") or "").strip()
        if not reason:
            violations.append(Violation("exclusion_reason", str(row["stable_id"]),
                                        "exclusion has no reason"))
        elif len(reason) < 10:
            violations.append(Violation("exclusion_reason", str(row["stable_id"]),
                                        f"reason too vague to act on: {reason!r}"))
    return violations


def check_glyph_pairing(resolution_rows: Iterable[dict]) -> List[Violation]:
    """A glyph id without its diagram is meaningless — ids are unique only
    within a diagram, so one without the other cannot be resolved."""
    violations = []
    for row in resolution_rows:
        glyph = str(row.get("glyph_id") or "").strip()
        diagram = str(row.get("diagram_stid") or "").strip()
        if bool(glyph) != bool(diagram):
            violations.append(Violation(
                "glyph_pairing", str(row["uuid"]),
                f"glyph_id={glyph!r} and diagram_stid={diagram!r} must be "
                f"present together or absent together"))
    return violations


def check_release_recorded(resolution_rows: Iterable[dict],
                           expected_release: str | None = None) -> List[Violation]:
    """Every row carries its Reactome release.

    Version skew has produced a false finding on this project; a mapping
    without a release cannot be refused when it does not match.
    """
    violations = []
    seen: Set[str] = set()
    for row in resolution_rows:
        release = str(row.get("release") or "").strip()
        if not release:
            violations.append(Violation("release_recorded", str(row["uuid"]),
                                        "row carries no Reactome release"))
        else:
            seen.add(release)
    if len(seen) > 1:
        violations.append(Violation("release_recorded", "<table>",
                                    f"mixed releases in one table: {sorted(seen)}"))
    if expected_release and seen and seen != {str(expected_release)}:
        violations.append(Violation("release_recorded", "<table>",
                                    f"expected release {expected_release}, found {sorted(seen)}"))
    return violations


def run_all(network_node_uuids: Set[str],
            resolution_rows: List[dict],
            exclusion_rows: List[dict],
            participating_entities: Set[str] | None = None,
            expected_release: str | None = None) -> List[Violation]:
    """All checks. Empty result means the mapping is complete in both
    directions; anything else names what is missing and why it matters."""
    violations: List[Violation] = []
    violations += check_reverse_completeness(network_node_uuids, resolution_rows)
    violations += check_referential_integrity(network_node_uuids, resolution_rows)
    violations += check_exclusions_have_reasons(exclusion_rows)
    violations += check_glyph_pairing(resolution_rows)
    violations += check_release_recorded(resolution_rows, expected_release)
    if participating_entities is not None:
        violations += check_forward_completeness(participating_entities,
                                                 resolution_rows, exclusion_rows)
    return violations


def summarise(violations: Iterable[Violation]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for v in violations:
        counts[v.check] = counts.get(v.check, 0) + 1
    return counts
