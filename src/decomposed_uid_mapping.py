import pandas as pd

decomposed_uid_mapping_column_types = {
    "uid": str,
    "reactome_id": int,
    "component_id": int,
    "component_id_or_reference_entity_id": pd.Int64Dtype(),
    "input_or_output_uid": str,
    "input_or_output_reactome_id": pd.Int64Dtype(),
    "stoichiometry": pd.Int64Dtype(),
    "entity_stoichiometry": pd.Int64Dtype(),  # For original entity stoichiometry
    "component_stoichiometry": pd.Int64Dtype(),  # For complex component stoichiometry
}
