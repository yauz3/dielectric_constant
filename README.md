# SMILES-Driven Dielectric Constant Prediction for Bioengineering-Oriented Solvent and Polymer Screening

This repository contains the **code and data** used in the study:

**SMILES-Driven Dielectric Constant Prediction for Bioengineering-Oriented Solvent and Polymer Screening**

The project focuses on building **SMILES-based machine learning pipelines** to predict the **dielectric constant (relative permittivity, εr)** for screening solvents and polymer-like chemistries in materials and bioengineering-oriented workflows.

---

## Contents

- Source code for data processing, model training, evaluation, and result reporting  
- Curated dataset(s) used in the experiments  
- Configuration files / scripts to reproduce the benchmark runs  
- (Optional) Trained models and prediction outputs (if included)

> If you use this repository, please cite the associated manuscript (see Citation section below).

---

## Project Structure

You may adapt the structure below to match your current files:

```text
.
├── data/
│   ├── raw/                   # original/raw data (if applicable)
│   ├── processed/             # cleaned / feature-ready data
│   └── README.md              # data notes (sources, preprocessing, licensing)
├── src/
│   ├── preprocessing/         # SMILES cleaning, featurization, scaling, etc.
│   ├── models/                # training scripts / model definitions
│   ├── evaluation/            # metrics, plots, reports
│   └── utils/                 # helper functions
├── notebooks/                 # exploratory analysis (optional)
├── results/                   # figures, tables, logs (optional)
├── requirements.txt
└── README.md
