"""Unit tests for the schema-backed provenance exports (nodes.csv,
node_reaction_context.csv). No Neo4j: entity labels / terminal components are
mocked. See schema/logic_network.linkml.yaml."""
import pandas as pd

import src.logic_network_generator as m
from src import neo4j_connector


def test_export_nodes_classifies_kinds(tmp_path, monkeypatch):
    # simple entity labels; variant set-derivation stubbed out (no Neo4j)
    monkeypatch.setattr(neo4j_connector, "get_labels",
                        lambda e: ["EntityWithAccessionedSequence"])
    monkeypatch.setattr(m, "_derive_sets_and_chosen",
                        lambda parent, members: (["R-HSA-75202"], ["R-HSA-68891"]))
    monkeypatch.setattr(m, "get_terminal_components", lambda s: {s})

    variant = "R-HSA-141608::variant::R-HSA-68365_R-HSA-68891"
    # Realistic UUIDs (export_nodes detects mapping direction by UUID shape).
    s1 = "aaaaaaaa-0000-0000-0000-000000000001"
    v1 = "aaaaaaaa-0000-0000-0000-000000000002"
    d1 = "aaaaaaaa-0000-0000-0000-000000000003"
    rxn1 = "aaaaaaaa-0000-0000-0000-0000000000r1"
    edges = pd.DataFrame([
        {"source_id": s1, "target_id": rxn1, "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": rxn1, "target_id": v1, "pos_neg": "pos", "and_or": None,
         "edge_type": "output", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": v1, "target_id": d1, "pos_neg": "pos", "and_or": "and",
         "edge_type": "dissociation", "stoichiometry": 1, "edge_reaction_id": None},
    ])
    reaction_id_map = pd.DataFrame({"uid": [rxn1], "reactome_id": ["R-HSA-100"]})
    uuid_mapping = {s1: "R-HSA-999", v1: variant, d1: "R-HSA-888"}

    out = tmp_path / "nodes.csv"
    m.export_nodes(edges, reaction_id_map, uuid_mapping, str(out))
    rows = {r["uuid"]: r for r in pd.read_csv(out, dtype=str, keep_default_na=False)
            .to_dict("records")}
    assert rows[rxn1]["node_kind"] == "reaction"
    assert rows[rxn1]["diagram_entity_id"] == "R-HSA-100"
    assert rows[v1]["node_kind"] == "set_variant"
    assert rows[v1]["diagram_entity_id"] == "R-HSA-141608"
    assert "R-HSA-68891" in rows[v1]["member_leaves"]
    assert rows[s1]["node_kind"] == "simple_entity"
    assert rows[d1]["node_kind"] == "dissociation_sink"


def test_export_node_reaction_context(tmp_path):
    """Catalyst rows name the wired member, not the fetch row's parent uuid.

    This test previously asserted the OPPOSITE — that the parent uuid "c1"
    from the catalyst fetch row appears in the export. That expectation is
    what kept #67 alive: the export was orphaning 100% of catalyst and
    regulator rows and a green test said it was correct. Kept here, inverted,
    as the record of what the right answer is.
    """
    entity_uuid_registry = {
        ("R-HSA-999", "rxn1", "input"): "s1",
        ("R-HSA-141608::variant::x", "rxn1", "output"): "v1",
    }
    reaction_id_map = pd.DataFrame({"uid": ["rxn1"], "reactome_id": ["R-HSA-100"]})
    # "c1" is the undecomposed parent the fetch row carries; "c1_member" is the
    # terminal member append_regulators actually wires to the reaction.
    catreg = pd.DataFrame([{"reaction_id": "R-HSA-100", "entity_id": "R-HSA-7",
                            "edge_type": "catalyst", "uuid": "c1", "reaction_uuid": "rxn1"}])
    edges = pd.DataFrame([
        {"source_id": "s1", "target_id": "rxn1", "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": "rxn1", "target_id": "v1", "pos_neg": "pos", "and_or": None,
         "edge_type": "output", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": "c1_member", "target_id": "rxn1", "pos_neg": "pos", "and_or": "and",
         "edge_type": "catalyst", "stoichiometry": 1, "edge_reaction_id": None},
    ])
    out = tmp_path / "ctx.csv"
    m.export_node_reaction_context(entity_uuid_registry, reaction_id_map, catreg,
                                   str(out), logic_network=edges)
    ctx = pd.read_csv(out).to_dict("records")
    triples = {(r["context_node"], r["reaction_id"], r["role"]) for r in ctx}
    assert ("s1", "R-HSA-100", "input") in triples
    assert ("v1", "R-HSA-100", "output") in triples
    assert ("c1_member", "R-HSA-100", "catalyst") in triples
    assert ("c1", "R-HSA-100", "catalyst") not in triples


def test_parse_variant_members():
    assert m._parse_variant_members("R-HSA-1::variant::R-HSA-2_R-HSA-3") == {"R-HSA-2", "R-HSA-3"}
    assert m._parse_variant_members("R-HSA-1") == set()


def test_context_export_emits_no_orphaned_rows(tmp_path):
    """Every context row must name a node that exists in the logic network.

    LNG #67: append_regulators DECOMPOSES a catalyst or regulator to its
    terminal members and emits edges from those member UUIDs, but the export
    wrote the UUID carried on the raw catalyst/regulator fetch row — the
    UNDECOMPOSED PARENT. Measured on the ten-pathway catalog before the fix:
    100% of 2,228 catalyst rows and 100% of 1,018 regulator rows named a node
    absent from logic_network.csv, 31.3% of the whole export, while input and
    output were clean at 0%.

    That silently broke every consumer that tries to answer "which node is
    this entity, at this reaction, in this role" for exactly the two roles
    that carry the causality.
    """
    rxn1 = "aaaaaaaa-0000-0000-0000-0000000000r1"
    inp = "aaaaaaaa-0000-0000-0000-000000000001"
    out_uuid = "aaaaaaaa-0000-0000-0000-000000000002"
    # The decomposed member the network actually wires up...
    member = "aaaaaaaa-0000-0000-0000-00000000000c"
    # ...versus the parent complex the fetch row carries. Distinct on purpose:
    # when they coincide the bug is invisible.
    parent = "aaaaaaaa-0000-0000-0000-0000000000ff"

    edges = pd.DataFrame([
        {"source_id": inp, "target_id": rxn1, "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": rxn1, "target_id": out_uuid, "pos_neg": "pos", "and_or": None,
         "edge_type": "output", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
        {"source_id": member, "target_id": rxn1, "pos_neg": "pos", "and_or": "and",
         "edge_type": "catalyst", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
    ])
    reaction_id_map = pd.DataFrame({"uid": [rxn1], "reactome_id": ["R-HSA-100"]})
    entity_uuid_registry = {
        ("R-HSA-11", rxn1, "input"): inp,
        ("R-HSA-22", rxn1, "output"): out_uuid,
    }
    catalyst_regulator_map = pd.DataFrame([
        {"uuid": parent, "reaction_id": "R-HSA-100", "edge_type": "catalyst"},
    ])

    out = tmp_path / "node_reaction_context.csv"
    m.export_node_reaction_context(entity_uuid_registry, reaction_id_map,
                                   catalyst_regulator_map, str(out),
                                   logic_network=edges)

    rows = pd.read_csv(out, dtype=str, keep_default_na=False).to_dict("records")
    live = set(edges["source_id"]) | set(edges["target_id"])
    orphans = [r for r in rows if r["context_node"] not in live]
    assert not orphans, f"orphaned context rows: {orphans}"

    # And the catalyst must be present at all — dropping the row instead of
    # fixing the UUID would also make the orphan count zero.
    catalysts = [r for r in rows if r["role"] == "catalyst"]
    assert len(catalysts) == 1, f"expected one catalyst row, got {catalysts}"
    assert catalysts[0]["context_node"] == member, (
        "catalyst row must name the decomposed member the network wires up, "
        f"not the parent fetch-row UUID; got {catalysts[0]['context_node']}"
    )


def test_glyph_id_and_diagram_are_written_together(tmp_path, monkeypatch):
    """A glyph id without its diagram is meaningless, and vice versa.

    Glyph ids are unique only WITHIN a diagram — the same integer identifies a
    different drawing in another one — so a row carrying one without the other
    cannot be resolved back to anything. This is the invariant behind the
    original question: knowing a uuid came from glyph 535 is only useful if you
    also know which diagram 535 belongs to.
    """
    monkeypatch.setattr(m, "get_labels", lambda e: ["EntityWithAccessionedSequence"], raising=False)

    rxn = "aaaaaaaa-0000-0000-0000-0000000000r1"
    src = "aaaaaaaa-0000-0000-0000-000000000001"
    edges = pd.DataFrame([
        {"source_id": src, "target_id": rxn, "pos_neg": "pos", "and_or": "and",
         "edge_type": "input", "stoichiometry": 1, "edge_reaction_id": "R-HSA-100"},
    ])
    reaction_id_map = pd.DataFrame({"uid": [rxn], "reactome_id": ["R-HSA-100"]})

    import src.diagram_connectivity as dc
    # One entity at one reaction, drawn once: the triple resolves to a glyph.
    monkeypatch.setattr(dc, "diagram_glyph_positions",
                        lambda pid: {("R-HSA-100", "R-HSA-999", "input"): [535]})
    monkeypatch.setattr(dc, "covering_diagram_stid", lambda pid: "R-HSA-1257604")
    monkeypatch.setattr(m, "get_pathway_participating_entities", lambda pid: set(), raising=False)

    out = tmp_path / "node_resolution.csv"
    exc = tmp_path / "node_exclusions.csv"
    m.export_node_resolution("R-HSA-100", edges, reaction_id_map,
                             {src: "R-HSA-999"}, str(out), str(exc))

    rows = pd.read_csv(out, dtype=str, keep_default_na=False).to_dict("records")
    for row in rows:
        has_glyph = bool(row["glyph_id"].strip())
        has_diagram = bool(row["diagram_stid"].strip())
        assert has_glyph == has_diagram, (
            f"glyph_id={row['glyph_id']!r} and diagram_stid={row['diagram_stid']!r} "
            "must be present together or absent together"
        )


def test_export_cofactors_lists_all_and_flags_present(tmp_path, monkeypatch):
    """Every known cofactor is listed; only those in the network are flagged.

    Listing all of them is what lets a consumer tell "this pathway has no
    cofactors" from "this bundle predates the file".
    """
    monkeypatch.setattr(neo4j_connector, "get_cofactor_species", lambda: [
        {"stable_id": "R-ALL-113592", "molecule": "ATP", "chebi_id": "30616",
         "name": "ATP [cytosol]"},
        {"stable_id": "R-ALL-29356", "molecule": "H2O", "chebi_id": "15377",
         "name": "H2O [cytosol]"},
    ])
    monkeypatch.setattr(neo4j_connector, "get_reactome_release", lambda: 97)

    edges = pd.DataFrame([{"source_id": "u-atp", "target_id": "u-rxn"}])
    out = tmp_path / "cofactors.csv"
    m.export_cofactors(edges, {"R-ALL-113592": "u-atp"}, str(out))

    df = pd.read_csv(out)
    assert list(df.columns) == ["stable_id", "molecule", "chebi_id", "name",
                                "in_network", "reactome_release"]
    assert len(df) == 2, "every known cofactor is listed, not only the present ones"
    assert set(df.loc[df.in_network == 1, "stable_id"]) == {"R-ALL-113592"}
    assert set(df["reactome_release"]) == {97}, "the release must travel with the list"


def test_export_cofactors_writes_a_file_even_when_none_are_present(tmp_path, monkeypatch):
    """A pathway with no cofactors still gets the file, all flags zero.

    A missing file and an empty intersection mean different things and a
    consumer must be able to distinguish them.
    """
    monkeypatch.setattr(neo4j_connector, "get_cofactor_species", lambda: [
        {"stable_id": "R-ALL-113592", "molecule": "ATP", "chebi_id": "30616",
         "name": "ATP [cytosol]"},
    ])
    monkeypatch.setattr(neo4j_connector, "get_reactome_release", lambda: 97)

    edges = pd.DataFrame([{"source_id": "u-x", "target_id": "u-y"}])
    out = tmp_path / "cofactors.csv"
    m.export_cofactors(edges, {"R-HSA-9999": "u-x"}, str(out))

    df = pd.read_csv(out)
    assert out.exists()
    assert len(df) == 1
    assert int(df.in_network.sum()) == 0


def test_export_cofactors_handles_both_mapping_directions(tmp_path, monkeypatch):
    """`reactome_id_to_uuid` is stored either direction depending on caller.

    The first version of this exporter assumed stable_id -> uuid. Given the
    other direction it marked EVERY row absent and shipped a file saying no
    pathway contains any cofactor. The original test constructed the mapping in
    the assumed direction, so it passed either way and could not catch this.
    """
    monkeypatch.setattr(neo4j_connector, "get_cofactor_species", lambda: [
        {"stable_id": "R-ALL-113592", "molecule": "ATP", "chebi_id": "30616",
         "name": "ATP [cytosol]"},
    ])
    monkeypatch.setattr(neo4j_connector, "get_reactome_release", lambda: 97)

    uuid = "aaaaaaaa-0000-0000-0000-000000000001"
    edges = pd.DataFrame([{"source_id": uuid, "target_id": "u-rxn"}])

    for direction, mapping in (
        ("stable_id -> uuid", {"R-ALL-113592": uuid}),
        ("uuid -> stable_id", {uuid: "R-ALL-113592"}),
    ):
        out = tmp_path / f"cofactors_{direction.split()[0]}.csv"
        m.export_cofactors(edges, mapping, str(out))
        df = pd.read_csv(out)
        assert int(df.in_network.sum()) == 1, f"ATP missed with {direction}"


def test_export_cofactors_finds_an_entity_split_across_uuids(tmp_path, monkeypatch):
    """One stable id routinely maps to several uuids (the silo).

    A dict keyed by stable id holds only one of them, so scanning the mapping
    by key undercounts a split entity. GPVI carries four separate GTP nodes.
    """
    monkeypatch.setattr(neo4j_connector, "get_cofactor_species", lambda: [
        {"stable_id": "R-ALL-29438", "molecule": "GTP", "chebi_id": "37565",
         "name": "GTP [cytosol]"},
    ])
    monkeypatch.setattr(neo4j_connector, "get_reactome_release", lambda: 97)

    u1 = "aaaaaaaa-0000-0000-0000-00000000000a"
    u2 = "aaaaaaaa-0000-0000-0000-00000000000b"
    # Only the SECOND occurrence appears in the network.
    edges = pd.DataFrame([{"source_id": u2, "target_id": "u-rxn"}])
    out = tmp_path / "cofactors.csv"
    m.export_cofactors(edges, {u1: "R-ALL-29438", u2: "R-ALL-29438"}, str(out))

    df = pd.read_csv(out)
    assert int(df.in_network.sum()) == 1, "split entity missed when only one uuid is used"
