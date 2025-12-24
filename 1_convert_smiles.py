#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu


import pandas as pd
import requests
import time
from urllib.parse import quote
import pubchempy as pcp

# ============================
#  Ayarlar
# ============================

BASE_SLEEP = 0.5          # temel bekleme (sn)
MAX_RETRIES = 5           # PubChem için maksimum tekrar
CACTUS_URL = "https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"


# ============================
#  Yardımcı fonksiyonlar
# ============================

def normalize_name(raw_name: str) -> str:
    """
    Kimyasal isimleri PubChem/Cactus için biraz normalize eder.
    Örn:
      '(+)-a-pinene' -> 'alpha-pinene'
      'α-pinene'     -> 'alpha-pinene'
    """
    if not isinstance(raw_name, str):
        return raw_name

    name = raw_name.strip()

    # Başındaki (+)- / (-)- gibi işaretleri temizle
    for prefix in ["(+)-", "(-)-", "(+)", "(-)"]:
        if name.startswith(prefix):
            name = name[len(prefix):].strip()

    # Yunanca harfleri metne çevir
    replacements = {
        "α": "alpha",
        "β": "beta",
        "γ": "gamma",
        "δ": "delta",
    }
    for old, new in replacements.items():
        name = name.replace(old, new)

    # 'a-pinene' gibi yazımları 'alpha-pinene' yap
    name = name.replace("a-pinene", "alpha-pinene")

    return name


# ============================
#  Resolver 1: PubChem (pubchempy) + retry/backoff
# ============================

def name_to_smiles_pubchem(name: str) -> str | None:
    """
    PubChem (pubchempy) kullanarak isimden Canonical SMILES çeker.
    503/ServerBusy durumunda exponential backoff ile tekrar dener.
    """
    if not isinstance(name, str) or not name.strip():
        return None

    norm_name = normalize_name(name)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            compounds = pcp.get_compounds(norm_name, 'name')
            if not compounds:
                print(f"[PubChem] Bulunamadı: {name!r} (normalized: {norm_name!r})")
                return None

            smiles = compounds[0].canonical_smiles
            print(f"[PubChem] {name!r} -> {smiles}")
            # Başarılıysa hafif bekleyip çık
            time.sleep(BASE_SLEEP)
            return smiles

        except Exception as e:
            msg = str(e)
            # 503 / ServerBusy durumunda tekrar dene
            if "ServerBusy" in msg or "503" in msg:
                wait = BASE_SLEEP * (2 ** (attempt - 1))  # 0.5,1,2,4,8,...
                print(f"[PubChem] ServerBusy/503 ({name!r}), attempt {attempt}/{MAX_RETRIES}, {wait:.1f}s bekleniyor...")
                time.sleep(wait)
                continue
            else:
                # Başka bir hata ise retry etmenin anlamı yok
                print(f"[PubChem] Hata ({name!r}): {e}")
                return None

    print(f"[PubChem] Maksimum deneme sayısı aşıldı, vazgeçildi: {name!r}")
    return None


# ============================
#  Resolver 2: NCI Cactus
# ============================

def name_to_smiles_cactus(name: str) -> str | None:
    """
    NCI Cactus servisinden isimden SMILES çeker.
    """
    if not isinstance(name, str) or not name.strip():
        return None

    norm_name = normalize_name(name)
    encoded_name = quote(norm_name)
    url = CACTUS_URL.format(name=encoded_name)

    try:
        r = requests.get(url, timeout=10)
        print(f"[Cactus] İstek: {url} -> status {r.status_code}")
        if r.status_code != 200:
            return None

        smiles = r.text.strip()
        if not smiles or "Not Found" in smiles:
            print(f"[Cactus] Bulunamadı: {name!r} (normalized: {norm_name!r})")
            return None

        time.sleep(BASE_SLEEP)
        print(f"[Cactus] {name!r} -> {smiles}")
        return smiles
    except Exception as e:
        print(f"[Cactus] Hata ({name!r}): {e}")
        return None


# ============================
#  Ana resolver: önce PubChem, olmazsa Cactus
# ============================

def name_to_smiles(name: str) -> str | None:
    """
    Önce PubChem (retry/backoff'lu), olmazsa Cactus dener.
    """
    smiles = name_to_smiles_pubchem(name)
    if smiles:
        return smiles

    smiles = name_to_smiles_cactus(name)
    if smiles:
        return smiles

    print(f"[WARN] SMILES bulunamadı: {name!r}")
    return None


# ============================
#  CSV işleme
# ============================

def add_smiles_to_csv(input_csv: str, output_csv: str):
    df = pd.read_csv(input_csv)

    if "Name" not in df.columns:
        raise ValueError("CSV'de 'Name' kolonu bulunamadı. Kolon adlarını kontrol et.")

    smiles_list = []
    total = len(df)

    for idx, row in df.iterrows():
        name = row["Name"]
        print(f"\n[{idx + 1}/{total}] İşleniyor: {name!r}")
        smiles = name_to_smiles(name)
        smiles_list.append(smiles)

    df["SMILES"] = smiles_list
    df.to_csv(output_csv, index=False)
    print(f"\nSMILES eklenmiş CSV kaydedildi: {output_csv}")


# ============================
#  Çalıştırma
# ============================

if __name__ == "__main__":
    input_csv_path = "Dielectric_constants.csv"
    output_csv_path = "dielectric_nd_with_smiles.csv"
    add_smiles_to_csv(input_csv_path, output_csv_path)

