import pandas as pd


decomposed_uid_mapping_column_types = {
    "uid": str,
    "reactome_id": int,
    "component_id": int,
    "component_id_or_reference_entity_id": pd.Int64Dtype(),  # Use pd.Int64Dtype() for nullable integers
    "input_or_output_uid": str,
    "input_or_output_reactome_id": pd.Int64Dtype(),  # Use pd.Int64Dtype() for nullable integers
}
