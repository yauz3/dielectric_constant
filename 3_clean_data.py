#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu

import pandas as pd

def clean_mixed_feature_columns(
    input_csv: str,
    output_csv: str,
    id_columns=None,
    sample_n: int = 5,
):
    if id_columns is None:
        id_columns = ["Name", "Formula", "Dielectric constant(εrel)", "nD", "SMILES"]

    df = pd.read_csv(input_csv)
    print("Tüm kolonlar:", df.columns.tolist())
    print("Kimlik (dokunma) kolonları:", id_columns)

    for col in df.columns:
        if col in id_columns:
            print(f"\n[SKIP] Kimlik kolonu atlandı: {col!r}")
            continue

        print(f"\n[CHECK] Kolon: {col!r}")

        # Orijinal kolonu al
        s = df[col]

        # Eğer tüm değerler NaN ise uğraşmaya gerek yok
        if s.notna().sum() == 0:
            print("  - Tüm değerler NaN, geçiliyor.")
            continue

        # Deneme amaçlı sayıya çevir (numeric_s: float + NaN)
        numeric_s = pd.to_numeric(s, errors="coerce")

        # En az bir sayıya çevrilebilmiş değer var mı?
        has_any_numeric = numeric_s.notna().any()

        # Orijinalde boş olmayanlardan hangileri NaN’a düşmüş?
        mask_original_notna = s.notna()
        mask_became_nan = mask_original_notna & numeric_s.isna()

        num_strings_turned_nan = mask_became_nan.sum()

        print(f"  - Orijinal non-NaN sayısı: {mask_original_notna.sum()}")
        print(f"  - Sayıya çevrilebilen (non-NaN) değer sayısı: {numeric_s.notna().sum()}")
        print(f"  - Sayıya çevrilemeyip NaN olan (muhtemel string/hata) değer sayısı: {num_strings_turned_nan}")

        # Karışık kolon kriteri:
        # 1) En az bir değer sayıya çevrilebilmeli
        # 2) En az bir değer çevrilemeyip NaN'a düşmeli (string/hata)
        if has_any_numeric and num_strings_turned_nan > 0:
            print(f"  [CLEAN] Karışık tip tespit edildi, kolonu numeric + NaN yapıyoruz: {col!r}")

            # Debug için string örnekleri göster
            string_examples = s[mask_became_nan].head(sample_n).tolist()
            numeric_examples = numeric_s[numeric_s.notna()].head(sample_n).tolist()
            print(f"    - Örnek string/hata değerler: {string_examples}")
            print(f"    - Örnek numeric değerler: {numeric_examples}")

            # Artık kolonu 'numeric_s' ile değiştiriyoruz
            df[col] = numeric_s
        else:
            print(f"  [KEEP] Bu kolon temiz kalıyor (tamamen string ya da tamamen numeric): {col!r}")

    df.to_csv(output_csv, index=False)
    print(f"\nTemizlenmiş dosya kaydedildi: {output_csv}")


if __name__ == "__main__":
    input_file = "dielectric_nd_with_smiles_with_feature.csv"
    output_file = "dielectric_nd_with_smiles_with_feature_cleaned.csv"

    clean_mixed_feature_columns(input_file, output_file)

