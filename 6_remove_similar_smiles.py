#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu



from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator  # <-- yeni API
import pandas as pd
import numpy as np

def smiles_to_mol(smiles):
    return Chem.MolFromSmiles(smiles) if pd.notna(smiles) else None

def get_morgan_generator(radius=2, nBits=2048, _cache={}):
    """
    Aynı radius/nBits için MorganGenerator tek kez oluşturulsun diye basit cache.
    """
    key = (radius, nBits)
    if key not in _cache:
        _cache[key] = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,
            fpSize=nBits
        )
    return _cache[key]

def mol_to_fp(mol, radius=2, nBits=2048):
    if mol is None:
        return None
    gen = get_morgan_generator(radius=radius, nBits=nBits)
    return gen.GetFingerprint(mol)

def diversity_filter(df, smiles_col="SMILES", tanimoto_thresh=0.9, radius=2, nBits=2048):
    """
    tanimoto_thresh üzerinde benzerlik gösteren moleküllerden 
    yalnızca birini tutar. (Greedy diversity selection)
    """
    mols = [smiles_to_mol(s) for s in df[smiles_col]]
    fps = [mol_to_fp(m, radius=radius, nBits=nBits) for m in mols]

    keep_indices = []
    kept_fps = []

    for idx, fp in enumerate(fps):
        if fp is None:
            continue

        too_similar = False
        for kfp in kept_fps:
            sim = DataStructs.TanimotoSimilarity(fp, kfp)
            if sim >= tanimoto_thresh:
                too_similar = True
                break

        if not too_similar:
            keep_indices.append(idx)
            kept_fps.append(fp)

    print(
        f"[Diversity filter] Kept {len(keep_indices)} / {len(df)} molecules "
        f"(threshold = {tanimoto_thresh})"
    )

    return df.iloc[keep_indices].reset_index(drop=True)

# Örnek kullanım:
df = pd.read_csv("dielectric_nd_with_smiles_with_feature_cleaned_2_with_fingerprint.csv")
df_filtered = diversity_filter(df, smiles_col="SMILES", tanimoto_thresh=0.90)
df_filtered.to_csv(
    "dielectric_nd_with_smiles_with_feature_cleaned_2_with_fingerprint_fingerprint_filtration.csv",
    index=False
)

