#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu

import pandas as pd

def drop_pure_string_columns(
    input_csv: str,
    output_csv: str,
    id_columns=None,
    sample_n: int = 5,
):
    if id_columns is None:
        id_columns = ["Name", "Formula", "Dielectric constant(εrel)", "nD", "SMILES"]

    df = pd.read_csv(input_csv)
    print("Tüm kolonlar:", df.columns.tolist())
    print("Kimlik (korunacak) kolonlar:", id_columns)

    cols_to_drop = []

    for col in df.columns:
        if col in id_columns:
            print(f"\n[SKIP] Kimlik kolonu, asla silme: {col!r}")
            continue

        print(f"\n[CHECK] Kolon: {col!r}")

        s = df[col]
        non_na = s.dropna()

        # Tamamen boş (tüm NaN) kolonlar:
        if non_na.empty:
            # İstersen burayı aktif edip tamamen boş kolonları da silebilirsin:
            # cols_to_drop.append(col)
            print("  - Bu kolon tamamen NaN (boş). Şimdilik dokunmuyoruz.")
            continue

        # Non-NaN değerlerin hepsi string mi?
        all_str = non_na.apply(lambda x: isinstance(x, str)).all()

        if all_str:
            print(f"  [DROP] Bu kolon sadece string (error/not) içeriyor, silinecek: {col!r}")
            print("    - Örnek string değerler:", non_na.head(sample_n).tolist())
            cols_to_drop.append(col)
        else:
            print(f"  [KEEP] Bu kolonda string dışında (muhtemelen sayısal) değerler de var, tutuluyor: {col!r}")
            print("    - Örnek non-NaN değerler:", non_na.head(sample_n).tolist())

    if cols_to_drop:
        print("\nSilinecek kolonlar:", cols_to_drop)
        df = df.drop(columns=cols_to_drop)
    else:
        print("\nSilinecek kolon bulunamadı (sadece string içeren yok).")

    df.to_csv(output_csv, index=False)
    print(f"\nFinal temizlenmiş dosya kaydedildi: {output_csv}")


if __name__ == "__main__":
    input_file = "dielectric_nd_with_smiles_with_feature_cleaned.csv"
    output_file = "dielectric_nd_with_smiles_with_feature_cleaned_2.csv"

    drop_pure_string_columns(input_file, output_file)

