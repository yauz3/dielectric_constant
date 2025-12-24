# SMILES-Driven Dielectric Constant Prediction for Bioengineering-Oriented Solvent and Polymer Screening

This repository contains the **code and datasets** used in the study:

**SMILES-Driven Dielectric Constant Prediction for Bioengineering-Oriented Solvent and Polymer Screening**

The project implements a **SMILES-based, leakage-aware machine learning workflow** for predicting the **dielectric constant (relative permittivity, εr)** to support **solvent and polymer-like chemistry screening** in materials and bioengineering-oriented contexts.

---

## Repository Contents

This repository includes:

- End-to-end preprocessing scripts (SMILES standardization, cleaning, deduplication)
- Fingerprint-based feature generation (RDKit required)
- Train/test split generation
- Model training and evaluation pipeline(s)
- Curated CSV datasets produced at different pipeline stages

---

## File Structure (as in this repository)

```text
.
├── 1_convert_smiles.py
├── 2_read_SMILES_to_features.py
├── 3_clean_data.py
├── 4_clean_data_str_coloums.py
├── 5_prepare_fingerprint_features.py
├── 6_remove_similar_smiles.py
├── 7_split_data_into_train_test.py
├── 9_model_train_and_validate.py
├── Dielectric_constants.csv
├── dielectric_nd_with_smiles_with_feature_cleaned_2_with_fingerprint_fingerprint_filtration.csv
├── train_ready.csv
└── test_ready.csv
