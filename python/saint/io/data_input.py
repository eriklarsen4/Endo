# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 23:36:08 2026

@author: Erik
"""

# %% Imports

import pandas as pd


# %% Wide to long converter with replicate collapse

def extract_bait_matrix(df):
    """
    Convert a wide format dataframe into the long format required by SAINT.

    The input df has a Protein column and additional columns of the form
    condition underscore bait underscore rep.

    Replicates for the same bait are collapsed by summing their values.

    The output is a long format dataframe with columns
    Bait, Protein, and Count_1.
    """

    protein = df["Protein"].tolist()
    value_cols = [c for c in df.columns if c != "Protein"]

    # Parse wide columns into (bait_name, numeric_values)
    parsed = {}
    for col in value_cols:
        parts = col.split("_")
        if len(parts) < 3:
            continue

        bait_name = parts[1]
        counts = pd.to_numeric(df[col], errors="coerce").fillna(0).to_numpy()

        if bait_name not in parsed:
            parsed[bait_name] = counts.copy()
        else:
            parsed[bait_name] += counts

    # Build long-format dataframe
    records = []
    for bait_name, collapsed_counts in parsed.items():
        for p, x in zip(protein, collapsed_counts):
            records.append(
                {
                    "Bait": bait_name,
                    "Protein": p,
                    "Count_1": float(x)
                }
            )

    long_df = pd.DataFrame(records)
    return long_df


