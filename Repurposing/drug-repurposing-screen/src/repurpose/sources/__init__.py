"""Data-source connectors.

Each connector returns a *canonical* pandas DataFrame so the pipeline does not
care whether the bytes came from the bundled demo CSVs, a ChEMBL SQLite dump,
Open Targets parquet, or a DrugBank XML export.

Canonical schemas
-----------------
drugs:            drug_id, drug_name, max_phase, approved_us, approved_eu, modality
drug_targets:     drug_id, target_symbol, action_type, mechanism_of_action
drug_indications: drug_id, efo_id, indication_name
target_disease:   target_symbol, efo_id, disease_name, assoc_score   (0..1)
disease_ontology: efo_id, disease_name, parent_efo_id
"""
