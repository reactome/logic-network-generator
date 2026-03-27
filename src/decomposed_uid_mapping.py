import pandas as pd

decomposed_uid_mapping_column_types = {
    "uid": str,
    "reactome_id": str,  # The reaction stId this entity participates in
    "component_id": str,
    "component_id_or_reference_entity_id": str,
    "input_or_output_uid": str,
    "input_or_output_reactome_id": str,
    "source_entity_id": str,  # The parent entity (Complex or EntitySet) that was decomposed
    "source_reaction_id": str,  # The original Reactome reaction (for virtual reactions)
    "stoichiometry": "Int64",  # Stoichiometric coefficient from Neo4j hasComponent relationships
}
