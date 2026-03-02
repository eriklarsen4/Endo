# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 23:36:08 2026

@author: Erik
"""

# %% Imports

import pandas as pd
from pathlib import Path

# %% 

def load_bait_data(input_obj):
    """
    Flexible loader that accepts either:
      - a pandas DataFrame
      - a file path to a CSV

    Returns:
      bait_list: list of bait names
      X_by_bait: dict mapping bait → numeric matrix
      metadata: dict containing bait-level biological metadata
    """

    # Case 1: user passed a DataFrame
    if isinstance(input_obj, pd.DataFrame):
        df = input_obj.copy()

    # Case 2: user passed a file path
    else:
        path = Path(input_obj)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {input_obj}")
        df = pd.read_csv(path)

    # Convert wide → long
    long_df = extract_bait_matrix(df)

    # Identify all baits
    bait_list = sorted(long_df["Bait"].unique())

    # Split into numeric matrices per bait
    X_by_bait = {}
    for bait in bait_list:
        df_b = long_df[long_df["Bait"] == bait].copy()

        # Sort by Protein to ensure consistent row order
        df_b = df_b.sort_values("Protein").reset_index(drop=True)

        # Extract all numeric replicate columns (exclude metadata)
        numeric_cols = [c for c in df_b.columns if c not in ("Protein", "Bait")]
        X_by_bait[bait] = df_b[numeric_cols].astype(float).to_numpy()

    metadata = {
    "baits": bait_list,
    "proteins_by_bait": {
        bait: df_b["Protein"].tolist()
        for bait, df_b in {
            bait: long_df[long_df["Bait"] == bait].sort_values("Protein")
            for bait in bait_list
        }.items()
    }
}

    return bait_list, X_by_bait, metadata

# %% Extract bait matrix (wide → long with replicate preservation)

def extract_bait_matrix(df):
    """
    Convert a wide format dataframe into the long format required by SAINT.

    The input df has a Protein column and additional columns of the form
    <condition>_<bait>_<rep>.

    Replicates are preserved as separate columns. No collapsing is performed
    in this function. Classical SAINT collapses replicates internally, while
    hierarchical SAINT uses the full replicate matrix.

    The output is a long format dataframe with columns:
    Protein, Bait, rep1, rep2, rep3, etc.
    """

    protein = df["Protein"].tolist()
    value_cols = [c for c in df.columns if c != "Protein"]

    # Parse wide columns into nested structure:
    # parsed[bait][rep_index] = vector of counts
    parsed = {}

    for col in value_cols:
        parts = col.split("_")
        if len(parts) < 3:
            continue

        bait_name = parts[1]
        rep_label = parts[2]

        try:
            rep_index = int(rep_label)
        except ValueError:
            continue

        counts = pd.to_numeric(df[col], errors="coerce").fillna(0).to_numpy()

        if bait_name not in parsed:
            parsed[bait_name] = {}

        parsed[bait_name][rep_index] = counts

    # Build long-format dataframe with replicate columns
    records = []

    for bait_name, rep_dict in parsed.items():
        max_rep = max(rep_dict.keys())

        for i, p in enumerate(protein):
            row = {
                "Protein": p,
                "Bait": bait_name
            }

            for r in range(1, max_rep + 1):
                row[f"rep{r}"] = float(rep_dict.get(r, pd.Series([0]*len(protein)))[i])

            records.append(row)

    long_df = pd.DataFrame(records)
    return long_df






