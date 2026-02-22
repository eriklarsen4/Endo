# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 23:36:08 2026

@author: Erik
"""

# %% Imports

import pandas as pd


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






