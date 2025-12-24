#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu


import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import ConvertToNumpyArray


def smiles_to_maccs(smiles_list):
    """
    Verilen SMILES listesi için MACCS fingerprint (167 bit) üretir
    ve bir DataFrame döner.
    """
    maccs_fps = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi) if pd.notna(smi) else None

        if mol is None:
            # Geçersiz SMILES için tümü 0 olan fingerprint
            maccs_arr = np.zeros((167,), dtype=int)
        else:
            fp = MACCSkeys.GenMACCSKeys(mol)
            maccs_arr = np.zeros((167,), dtype=int)
            ConvertToNumpyArray(fp, maccs_arr)

        maccs_fps.append(maccs_arr)

    df_maccs = pd.DataFrame(
        maccs_fps,
        columns=[f"MACCS_FP_{i}" for i in range(167)]
    )

    return df_maccs


def main():
    # 📌 Orijinal veri setini oku
    input_path = "dielectric_nd_with_smiles_with_feature_cleaned_2.csv"
    output_path = "dielectric_nd_with_smiles_with_feature_cleaned_2_with_fingerprint.csv"

    df = pd.read_csv(input_path)

    if "SMILES" not in df.columns:
        raise ValueError("Veri setinde 'SMILES' kolonu bulunamadı!")

    # 📌 SMILES'ları al
    smiles_list = df["SMILES"]

    # 📌 MACCS fingerprint üret
    df_maccs = smiles_to_maccs(smiles_list)

    # 📌 Orijinal veriyle birleştir
    df_final = pd.concat(
        [df.reset_index(drop=True), df_maccs.reset_index(drop=True)],
        axis=1
    )

    # 📌 Kaydet
    df_final.to_csv(output_path, index=False)

    print("✅ MACCS fingerprintleri eklendi ve kaydedildi:")
    print(f"   Girdi dosyası : {input_path}")
    print(f"   Çıktı dosyası : {output_path}")
    print(f"   Toplam satır  : {len(df_final)}")
    print(f"   Toplam kolon  : {df_final.shape[1]}")


if __name__ == "__main__":
    main()

