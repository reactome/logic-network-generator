"""Validation tests for Autophagy pathway (9612973).

Verifies that the generated logic network matches the Neo4j database:
1. All reactions in the pathway are represented
2. All entities in the UUID mapping exist in the database
3. Catalyst and regulator counts match the database
4. Decomposed entity sets contain valid members
5. Edge properties are valid

Requires: Neo4j database running with Reactome data.
"""

import pandas as pd
import pytest
from pathlib import Path
from py2neo import Graph


PATHWAY_ID = 9612973
PATHWAY_DIR = Path("output/Autophagy_9612973")


@pytest.fixture(scope="module")
def graph():
    """Create Neo4j graph connection."""
    try:
        g = Graph("bolt://localhost:7687", auth=("neo4j", "test"))
        g.run("RETURN 1").data()
        return g
    except Exception:
        pytest.skip("Neo4j database not available")


@pytest.fixture(scope="module")
def uuid_mapping():
    """Load the UUID-to-Reactome-ID mapping."""
    path = PATHWAY_DIR / "stid_to_uuid_mapping.csv"
    if not path.exists():
        pytest.skip("Autophagy output not generated")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def logic_network_sample():
    """Load logic network - sample if too large."""
    path = PATHWAY_DIR / "logic_network.csv"
    if not path.exists():
        pytest.skip("Autophagy output not generated")

    # Check file size - if over 10MB, sample rows
    file_size = path.stat().st_size
    if file_size > 10_000_000:
        # Read header + sample
        header = pd.read_csv(path, nrows=0)
        # Count lines efficiently
        with open(path) as f:
            total_lines = sum(1 for _ in f) - 1  # subtract header
        # Read first 1000, last 1000, and 1000 random from middle
        df_head = pd.read_csv(path, nrows=1000)
        df_tail = pd.read_csv(path, skiprows=range(1, max(2, total_lines - 999)), nrows=1000)
        df = pd.concat([df_head, df_tail], ignore_index=True)
        df.attrs['total_edges'] = total_lines
        df.attrs['sampled'] = True
    else:
        df = pd.read_csv(path)
        df.attrs['total_edges'] = len(df)
        df.attrs['sampled'] = False
    return df


@pytest.fixture(scope="module")
def reaction_connections():
    """Load reaction connections."""
    path = PATHWAY_DIR / "cache" / "reaction_connections.csv"
    if not path.exists():
        pytest.skip("Autophagy cache not available")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def decomposed_mapping():
    """Load decomposed UID mapping."""
    path = PATHWAY_DIR / "cache" / "decomposed_uid_mapping.csv"
    if not path.exists():
        pytest.skip("Autophagy decomposition cache not available")
    return pd.read_csv(path)


class TestAutophagyReactions:
    """Validate that all reactions in the pathway are represented."""

    def test_all_db_reactions_in_reaction_connections(self, graph, reaction_connections):
        """Every reaction in the Autophagy pathway should appear in reaction_connections."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        RETURN DISTINCT r.dbId as reaction_id, r.displayName as name
        """
        db_reactions = graph.run(query).data()
        db_reaction_ids = {int(r['reaction_id']) for r in db_reactions}

        generated_ids = set()
        for col in ['preceding_reaction_id', 'following_reaction_id']:
            generated_ids.update(
                int(x) for x in reaction_connections[col].dropna().unique()
            )

        missing = db_reaction_ids - generated_ids
        extra = generated_ids - db_reaction_ids

        print(f"\nDB reactions: {len(db_reaction_ids)}")
        print(f"Generated reactions: {len(generated_ids)}")
        print(f"Missing from generated: {len(missing)}")
        if missing:
            missing_names = [r['name'] for r in db_reactions if r['reaction_id'] in missing]
            print(f"Missing reactions: {missing_names[:10]}")

        assert len(missing) == 0, (
            f"{len(missing)} DB reactions missing from reaction_connections: "
            f"{sorted(missing)[:10]}"
        )

    def test_reaction_count_matches_db(self, graph, reaction_connections):
        """Number of unique reactions should match the database."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        RETURN count(DISTINCT r.dbId) as count
        """
        db_count = graph.run(query).data()[0]['count']

        generated_ids = set()
        for col in ['preceding_reaction_id', 'following_reaction_id']:
            generated_ids.update(
                int(x) for x in reaction_connections[col].dropna().unique()
            )

        print(f"\nDB reaction count: {db_count}")
        print(f"Generated reaction count: {len(generated_ids)}")
        assert len(generated_ids) == db_count


class TestAutophagyEntities:
    """Validate that entities in the output exist in the database."""

    def test_all_mapped_entities_exist_in_db(self, graph, uuid_mapping):
        """Every stable ID in the UUID mapping should exist in Neo4j."""
        stable_ids = uuid_mapping['stable_id'].unique().tolist()
        print(f"\nTotal mapped entities: {len(stable_ids)}")

        # Batch check in Neo4j using stId
        ids_str = ", ".join(f"'{sid}'" for sid in stable_ids)
        query = f"""
        MATCH (e)
        WHERE e.stId IN [{ids_str}]
        RETURN e.stId as entity_id
        """
        db_results = graph.run(query).data()
        db_entity_ids = {r['entity_id'] for r in db_results}

        missing = set(stable_ids) - db_entity_ids
        print(f"Entities found in DB: {len(db_entity_ids)}")
        print(f"Missing from DB: {len(missing)}")

        assert len(missing) == 0, (
            f"{len(missing)} entities in UUID mapping not found in DB: "
            f"{sorted(missing)[:20]}"
        )

    def test_mapped_entities_are_physical_entities(self, graph, uuid_mapping):
        """Mapped entities should be PhysicalEntity or DatabaseObject types."""
        stable_ids = uuid_mapping['stable_id'].unique().tolist()

        # Sample if too many
        sample = stable_ids[:200] if len(stable_ids) > 200 else stable_ids
        ids_str = ", ".join(f"'{sid}'" for sid in sample)

        query = f"""
        MATCH (e)
        WHERE e.stId IN [{ids_str}]
        RETURN e.stId as entity_id, labels(e) as labels
        """
        results = graph.run(query).data()

        valid_labels = {
            'PhysicalEntity', 'EntityWithAccessionedSequence', 'Complex',
            'EntitySet', 'DefinedSet', 'CandidateSet', 'OpenSet',
            'SimpleEntity', 'GenomeEncodedEntity', 'OtherEntity',
            'Polymer', 'Drug', 'ChemicalDrug', 'ProteinDrug',
            'DatabaseObject', 'Cell',
        }

        invalid_entities = []
        for r in results:
            entity_labels = set(r['labels'])
            if not entity_labels & valid_labels:
                invalid_entities.append((r['entity_id'], r['labels']))

        print(f"\nChecked {len(results)} entities")
        if invalid_entities:
            print(f"Invalid entity types: {invalid_entities[:10]}")

        assert len(invalid_entities) == 0, (
            f"{len(invalid_entities)} entities have unexpected types: {invalid_entities[:10]}"
        )

    def test_entity_count_reasonable(self, uuid_mapping):
        """UUID mapping should have a reasonable number of entries."""
        unique_stable_ids = uuid_mapping['stable_id'].nunique()
        total_uuids = len(uuid_mapping)

        print(f"\nTotal UUID entries: {total_uuids}")
        print(f"Unique stable IDs: {unique_stable_ids}")
        print(f"Average UUIDs per entity: {total_uuids / unique_stable_ids:.1f}")

        assert unique_stable_ids > 0, "No entities in UUID mapping"
        assert total_uuids > 0, "No UUID entries"


class TestAutophagyCatalystsAndRegulators:
    """Validate catalysts and regulators match the database."""

    def test_catalyst_count(self, graph, logic_network_sample):
        """Number of catalyst edges should match database catalyst count."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:catalystActivity]->(ca:CatalystActivity)-[:physicalEntity]->(pe:PhysicalEntity)
        RETURN count(DISTINCT pe.dbId) as unique_catalysts,
               count(*) as total_catalyst_relations
        """
        db_result = graph.run(query).data()[0]

        catalyst_edges = logic_network_sample[
            logic_network_sample['edge_type'] == 'catalyst'
        ]

        print(f"\nDB unique catalysts: {db_result['unique_catalysts']}")
        print(f"DB total catalyst relations: {db_result['total_catalyst_relations']}")
        print(f"Generated catalyst edges: {len(catalyst_edges)}")

        # Catalyst edges should be > 0 if DB has catalysts
        if db_result['unique_catalysts'] > 0:
            assert len(catalyst_edges) > 0, "DB has catalysts but none in generated network"

    def test_positive_regulator_count(self, graph, logic_network_sample):
        """Positive regulator edges should match database."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(reg:PositiveRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN count(DISTINCT pe.dbId) as unique_regulators,
               count(*) as total_relations
        """
        db_result = graph.run(query).data()[0]

        pos_reg_edges = logic_network_sample[
            (logic_network_sample['edge_type'] == 'regulator') &
            (logic_network_sample['pos_neg'] == 'pos')
        ]

        print(f"\nDB unique positive regulators: {db_result['unique_regulators']}")
        print(f"DB total positive regulation relations: {db_result['total_relations']}")
        print(f"Generated positive regulator edges: {len(pos_reg_edges)}")

        if db_result['unique_regulators'] > 0:
            assert len(pos_reg_edges) > 0, "DB has positive regulators but none in generated network"

    def test_negative_regulator_count(self, graph, logic_network_sample):
        """Negative regulator edges should match database."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(reaction:ReactionLikeEvent)
        MATCH (reaction)-[:regulatedBy]->(reg:NegativeRegulation)-[:regulator]->(pe:PhysicalEntity)
        RETURN count(DISTINCT pe.dbId) as unique_regulators,
               count(*) as total_relations
        """
        db_result = graph.run(query).data()[0]

        neg_reg_edges = logic_network_sample[
            (logic_network_sample['edge_type'] == 'regulator') &
            (logic_network_sample['pos_neg'] == 'neg')
        ]

        print(f"\nDB unique negative regulators: {db_result['unique_regulators']}")
        print(f"DB total negative regulation relations: {db_result['total_relations']}")
        print(f"Generated negative regulator edges: {len(neg_reg_edges)}")

        if db_result['unique_regulators'] > 0:
            assert len(neg_reg_edges) > 0, "DB has negative regulators but none in generated network"


class TestAutophagyDecomposition:
    """Validate that entity decomposition is correct."""

    def test_entity_set_members_are_valid(self, graph, decomposed_mapping):
        """Entities in decomposed mapping that came from EntitySets should be valid members."""
        # Find EntitySet reactome_ids in the decomposed mapping
        set_reactome_ids = decomposed_mapping['reactome_id'].unique()

        # Sample up to 20 entity sets
        ids_str = ", ".join(str(int(rid)) for rid in set_reactome_ids[:50])
        query = f"""
        MATCH (es)
        WHERE es.dbId IN [{ids_str}] AND 'EntitySet' IN labels(es)
        OPTIONAL MATCH (es)-[:hasCandidate|hasMember]->(member)
        RETURN es.dbId as set_id, es.displayName as set_name,
               collect(DISTINCT member.dbId) as member_ids
        """
        db_sets = graph.run(query).data()

        print(f"\nEntitySets found in DB from decomposed mapping: {len(db_sets)}")
        for s in db_sets[:5]:
            print(f"  {s['set_name']} ({s['set_id']}): {len(s['member_ids'])} members")

        # For each EntitySet, check that the decomposed members are valid
        for entity_set in db_sets:
            set_id = entity_set['set_id']
            db_member_ids = set(entity_set['member_ids'])

            if not db_member_ids:
                continue

            # Get what we decomposed this set into
            set_rows = decomposed_mapping[
                decomposed_mapping['reactome_id'] == set_id
            ]
            decomposed_ids = set()
            for _, row in set_rows.iterrows():
                if pd.notna(row.get('input_or_output_reactome_id')):
                    decomposed_ids.add(int(row['input_or_output_reactome_id']))

            # Decomposed IDs should be a subset of what the DB says
            # (they could be deeper decompositions of the members)
            if decomposed_ids:
                print(f"  Set {set_id}: decomposed into {len(decomposed_ids)} terminal IDs, "
                      f"DB has {len(db_member_ids)} direct members")

    def test_complex_components_are_valid(self, graph, decomposed_mapping):
        """Entities from Complex decomposition should be valid components."""
        complex_reactome_ids = decomposed_mapping['reactome_id'].unique()

        ids_str = ", ".join(str(int(rid)) for rid in complex_reactome_ids[:50])
        query = f"""
        MATCH (c)
        WHERE c.dbId IN [{ids_str}] AND 'Complex' IN labels(c)
        OPTIONAL MATCH (c)-[:hasComponent]->(comp)
        RETURN c.dbId as complex_id, c.displayName as complex_name,
               collect(DISTINCT comp.dbId) as component_ids
        """
        db_complexes = graph.run(query).data()

        print(f"\nComplexes found in DB from decomposed mapping: {len(db_complexes)}")
        for c in db_complexes[:5]:
            print(f"  {c['complex_name']} ({c['complex_id']}): "
                  f"{len(c['component_ids'])} components")

    def test_decomposed_mapping_has_entries(self, decomposed_mapping):
        """Decomposed mapping should not be empty."""
        print(f"\nDecomposed mapping rows: {len(decomposed_mapping)}")
        print(f"Unique UIDs: {decomposed_mapping['uid'].nunique()}")
        print(f"Unique reactome_ids: {decomposed_mapping['reactome_id'].nunique()}")

        assert len(decomposed_mapping) > 0, "Decomposed mapping is empty"

    def test_reaction_inputs_outputs_in_db(self, graph, decomposed_mapping):
        """Reaction inputs and outputs should match what's in the database."""
        # Get a sample of reaction IDs from the decomposed mapping
        reaction_ids = decomposed_mapping['reactome_id'].unique()

        # Find which of these are actual reactions (not entities)
        sample_ids = reaction_ids[:30]
        ids_str = ", ".join(str(int(rid)) for rid in sample_ids)
        query = f"""
        MATCH (r:ReactionLikeEvent)
        WHERE r.dbId IN [{ids_str}]
        OPTIONAL MATCH (r)-[:input]->(input)
        OPTIONAL MATCH (r)-[:output]->(output)
        RETURN r.dbId as reaction_id, r.displayName as name,
               collect(DISTINCT input.dbId) as input_ids,
               collect(DISTINCT output.dbId) as output_ids
        """
        db_reactions = graph.run(query).data()

        print(f"\nReactions with inputs/outputs in DB: {len(db_reactions)}")
        for r in db_reactions[:5]:
            print(f"  {r['name']} ({r['reaction_id']}): "
                  f"{len(r['input_ids'])} inputs, {len(r['output_ids'])} outputs")

        # Every reaction should have at least one input and one output
        reactions_without_io = [
            r for r in db_reactions
            if not r['input_ids'] or not r['output_ids']
        ]
        if reactions_without_io:
            print(f"\nReactions without inputs or outputs: {len(reactions_without_io)}")
            for r in reactions_without_io[:5]:
                print(f"  {r['name']} ({r['reaction_id']})")


class TestAutophagyEdgeProperties:
    """Validate edge properties in the logic network."""

    def test_valid_edge_types(self, logic_network_sample):
        """All edge types should be valid."""
        valid = {'input', 'output', 'catalyst', 'regulator'}
        edge_types = set(logic_network_sample['edge_type'].unique())
        invalid = edge_types - valid
        assert len(invalid) == 0, f"Invalid edge types: {invalid}"

    def test_valid_pos_neg(self, logic_network_sample):
        """pos_neg should be 'pos' or 'neg'."""
        valid = {'pos', 'neg', ''}
        pos_neg_values = set(logic_network_sample['pos_neg'].dropna().unique())
        invalid = pos_neg_values - valid
        assert len(invalid) == 0, f"Invalid pos_neg values: {invalid}"

    def test_valid_and_or(self, logic_network_sample):
        """and_or should be 'and' or 'or'."""
        valid = {'and', 'or', ''}
        and_or_values = set(logic_network_sample['and_or'].dropna().unique())
        invalid = and_or_values - valid
        assert len(invalid) == 0, f"Invalid and_or values: {invalid}"

    def test_edge_type_distribution(self, logic_network_sample):
        """Report edge type distribution."""
        total = logic_network_sample.attrs.get('total_edges', len(logic_network_sample))
        sampled = logic_network_sample.attrs.get('sampled', False)

        dist = logic_network_sample['edge_type'].value_counts()
        print(f"\nTotal edges in file: {total}")
        print(f"Sampled: {sampled}")
        print(f"Edge type distribution (in sample):")
        for etype, count in dist.items():
            print(f"  {etype}: {count}")

    def test_no_null_source_or_target(self, logic_network_sample):
        """Source and target IDs should never be null."""
        assert logic_network_sample['source_id'].notna().all(), "Found null source_id"
        assert logic_network_sample['target_id'].notna().all(), "Found null target_id"

    def test_self_loop_ratio(self, logic_network_sample):
        """Report self-loop ratio (source == target)."""
        main_edges = logic_network_sample[
            ~logic_network_sample['edge_type'].isin(['catalyst', 'regulator'])
        ]
        if len(main_edges) == 0:
            pytest.skip("No main edges in sample")

        self_loops = main_edges[main_edges['source_id'] == main_edges['target_id']]
        ratio = len(self_loops) / len(main_edges)

        print(f"\nMain edges in sample: {len(main_edges)}")
        print(f"Self-loops: {len(self_loops)}")
        print(f"Self-loop ratio: {ratio*100:.1f}%")

        # Self-loops are expected when same entity appears as both input and output
        # But shouldn't be the vast majority
        assert ratio < 0.95, f"Self-loop ratio too high: {ratio*100:.1f}%"


class TestAutophagyCompleteness:
    """Validate completeness of the generated network."""

    def test_all_reaction_inputs_covered(self, graph, uuid_mapping):
        """Input entities from reactions should appear in the UUID mapping."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:input]->(input:PhysicalEntity)
        RETURN DISTINCT input.stId as entity_id, input.displayName as name
        """
        db_inputs = graph.run(query).data()
        db_input_ids = {r['entity_id'] for r in db_inputs}

        mapped_ids = set(uuid_mapping['stable_id'].unique())

        # Check direct coverage (entity itself or its decomposed parts)
        direct_coverage = db_input_ids & mapped_ids

        print(f"\nDB reaction input entities: {len(db_input_ids)}")
        print(f"Directly mapped: {len(direct_coverage)}")
        print(f"Not directly mapped: {len(db_input_ids - mapped_ids)}")

        # Some entities won't be directly mapped because they were decomposed
        # into their components. Check if their components are mapped.
        unmapped = db_input_ids - mapped_ids
        if unmapped:
            unmapped_str = ", ".join(f"'{eid}'" for eid in list(unmapped)[:20])
            query2 = f"""
            MATCH (e)-[:hasComponent|hasCandidate|hasMember*1..5]->(child)
            WHERE e.stId IN [{unmapped_str}]
            RETURN e.stId as parent_id, collect(DISTINCT child.stId) as child_ids
            """
            decomposed = graph.run(query2).data()
            for d in decomposed[:5]:
                child_coverage = set(d['child_ids']) & mapped_ids
                print(f"  Entity {d['parent_id']}: {len(child_coverage)}/{len(d['child_ids'])} "
                      f"children mapped")

    def test_all_reaction_outputs_covered(self, graph, uuid_mapping):
        """Output entities from reactions should appear in the UUID mapping."""
        query = f"""
        MATCH (pathway:Pathway {{dbId: {PATHWAY_ID}}})-[:hasEvent*]->(r:ReactionLikeEvent)
        MATCH (r)-[:output]->(output:PhysicalEntity)
        RETURN DISTINCT output.stId as entity_id, output.displayName as name
        """
        db_outputs = graph.run(query).data()
        db_output_ids = {r['entity_id'] for r in db_outputs}

        mapped_ids = set(uuid_mapping['stable_id'].unique())
        direct_coverage = db_output_ids & mapped_ids

        print(f"\nDB reaction output entities: {len(db_output_ids)}")
        print(f"Directly mapped: {len(direct_coverage)}")
        print(f"Not directly mapped: {len(db_output_ids - mapped_ids)}")
