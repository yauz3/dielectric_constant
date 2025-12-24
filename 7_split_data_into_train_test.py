#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu



import pandas as pd
from sklearn.model_selection import train_test_split

# 📌 Girdi dosyası
input_path = "dielectric_nd_with_smiles_with_feature_cleaned_2_with_fingerprint.csv"

# 📌 Çıktı dosyaları
train_path = "train_ready.csv"
test_path = "test_ready.csv"

# 1) Veriyi oku
df = pd.read_csv(input_path)
print(f"Toplam satır sayısı: {len(df)}")

# 2) Train / Test olarak ayır (80% / 20%)
train_df, test_df = train_test_split(
    df,
    test_size=0.2,
    random_state=42,
    shuffle=True
)

print(f"Train satır sayısı: {len(train_df)}")
print(f"Test satır sayısı : {len(test_df)}")

# 3) Indexleri sıfırla (opsiyonel ama temiz olur)
train_df = train_df.reset_index(drop=True)
test_df = test_df.reset_index(drop=True)

# 4) Kaydet
train_df.to_csv(train_path, index=False)
test_df.to_csv(test_path, index=False)

print(f"\n✅ Kaydedildi:")
print(f"  - {train_path}")
print(f"  - {test_path}")
