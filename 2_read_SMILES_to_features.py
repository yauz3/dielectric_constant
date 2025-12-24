#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/12/2025
# Author: Sadettin Y. Ugurlu

import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from mordred import Calculator, descriptors
import PyBioMed
from PyBioMed.PyMolecule import moe
from jpype import startJVM, shutdownJVM, JClass
import jpype
import os

print("The import has been successfully finished!")

current_dir = os.path.dirname(os.path.abspath(__file__))

# Set the path to the CDK JAR file
CDK_JAR_PATH = os.path.expanduser(f"{current_dir}/bin/cdk.jar")
print("CDK JAR path:", CDK_JAR_PATH)

# Start the JVM with CDK JAR
if not jpype.isJVMStarted():
    startJVM(jpype.getDefaultJVMPath(), "-ea", f"-Djava.class.path={CDK_JAR_PATH}")

# Load CDK classes
CDKDescCalc = JClass("org.openscience.cdk.qsar.DescriptorEngine")

# Load Mordred Calculator
mordred_calc = Calculator(descriptors, ignore_3D=True)


# 📌 Step 1: Read SMILES from CSV
def read_smiles(csv_file, smiles_column="SMILES"):
    """
    CSV dosyasını okur, SMILES kolonundan RDKit Mol objesi üretir.
    """
    df = pd.read_csv(csv_file)
    if smiles_column not in df.columns:
        raise ValueError(f"CSV'de '{smiles_column}' kolonu bulunamadı. Kolonları kontrol et: {df.columns.tolist()}")

    df["Mol"] = df[smiles_column].apply(
        lambda x: Chem.MolFromSmiles(x) if pd.notna(x) and isinstance(x, str) and x.strip() else None
    )
    return df


# 📌 Step 2: Compute RDKit Descriptors
def compute_rdkit_features(mol):
    if mol is None:
        return {}

    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
    }


# 📌 Step 3: Compute Mordred Descriptors
def compute_mordred_features(mol):
    if mol is None:
        return {}

    try:
        mordred_features = mordred_calc(mol)
        return dict(mordred_features)
    except Exception:
        return {}


# 📌 Step 4: Compute PyBioMed Descriptors
def compute_pybiomed_features(mol):
    """
    PyBioMed.MOE descriptor'larını hesaplar.
    Burada 'mol' bir RDKit Mol objesi olarak beklenir.
    """
    if mol is None:
        return {}

    try:
        mol_des = moe.GetMOE(mol)
        return mol_des
    except Exception:
        return {}


# 📌 Step 5: Compute CDK Descriptors
def compute_cdk_features(smiles):
    if smiles is None or not isinstance(smiles, str) or not smiles.strip():
        return {}

    try:
        # CDK DescriptorEngine, descriptor sınıfları listesi alır; burada örnek olarak sadece MW kullanıyoruz.
        cdk_calculator = CDKDescCalc("org.openscience.cdk.qsar.descriptors.molecular.MolecularWeightDescriptor")
        # DescriptorEngine normalde IAtomContainer listesi ile çalışır.
        # Basit kullanım için SMILES string'i üzerinden bir test yapılmış durumda;
        # ileride gerekirse buraya tam CDK IAtomContainer yaratımı eklenebilir.
        cdk_descriptors = cdk_calculator.calculate(smiles)
        return {"CDK_MolecularWeight": cdk_descriptors}
    except Exception:
        return {}


# 📌 Step 6: Process Dataset & Save Features
def generate_features(csv_file, output_csv, smiles_column="SMILES"):
    """
    dielectric_nd_with_smiles.csv dosyasını okuyup
    RDKit + Mordred + PyBioMed + CDK feature'larını hesaplar
    ve dielectric_nd_with_smiles_with_feature.csv olarak kaydeder.
    """
    df = read_smiles(csv_file, smiles_column=smiles_column)

    # Compute all descriptors
    print("RDKit descriptors are being calculated...")
    df["RDKit_Features"] = df["Mol"].apply(compute_rdkit_features)

    print("Mordred descriptors are being calculated...")
    df["Mordred_Features"] = df["Mol"].apply(compute_mordred_features)

    print("PyBioMed MOE descriptors are being calculated...")
    df["PyBioMed_Features"] = df["Mol"].apply(compute_pybiomed_features)

    print("CDK descriptors are being calculated...")
    df["CDK_Features"] = df[smiles_column].apply(compute_cdk_features)

    # Flatten feature dictionaries
    rdkit_df = pd.DataFrame(df["RDKit_Features"].tolist())
    mordred_df = pd.DataFrame(df["Mordred_Features"].tolist())
    pybiomed_df = pd.DataFrame(df["PyBioMed_Features"].tolist())
    cdk_df = pd.DataFrame(df["CDK_Features"].tolist())

    # Ana bilgileri de koruyalım (örn. Name, Formula, dielectric, nD, SMILES)
    base_cols = [c for c in df.columns if c not in ["Mol", "RDKit_Features", "Mordred_Features",
                                                    "PyBioMed_Features", "CDK_Features"]]
    base_df = df[base_cols].reset_index(drop=True)

    # Merge all features into final dataset
    final_df = pd.concat([base_df, rdkit_df, mordred_df, pybiomed_df, cdk_df], axis=1)

    # Save results
    final_df.to_csv(output_csv, index=False)
    print(f"Feature extraction complete. File saved as {output_csv}")


def merge_csv_on_smiles(file1, file2, output_file="merged_output.csv", smiles_column="SMILES"):
    """
    İki CSV dosyasını 'SMILES' kolonu üzerinden birleştirir (left-join, file1 ana dataset).
    """
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    for i, row in df2.iterrows():
        for col in df2.columns:
            try:
                df2.at[i, col] = float(df2.at[i, col])
            except Exception:
                if isinstance(row[col], str) and col != smiles_column:
                    df2.at[i, col] = None

    merged_df = df1.merge(df2, on=smiles_column, how="left")
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file saved as {output_file}")
    return merged_df


if __name__ == "__main__":
    input_csv_path = "dielectric_nd_with_smiles.csv"
    output_csv_path = "dielectric_nd_with_smiles_with_feature.csv"

    generate_features(input_csv_path, output_csv_path, smiles_column="SMILES")

    # İş bittiyse JVM'i kapatmak istersen:
    # if jpype.isJVMStarted():
    #     shutdownJVM()

