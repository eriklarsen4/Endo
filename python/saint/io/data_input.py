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
    Load IPMS data and construct per-bait replicate matrices for the hierarchical
    SAINT pipeline. This function accepts either a wide-format DataFrame or a
    path to a CSV file, converts the data into the long format required by the
    model, and assembles the per-bait numeric matrices used by the hierarchical
    EM algorithm and tau grid search.

    The loader preserves all replicates, performs no collapsing, and enforces a
    deterministic protein ordering within each bait. The resulting metadata
    includes the bait list and the per-bait protein lists, which define the
    canonical row ordering for all downstream model outputs.

    Parameters
    
    input_obj : DataFrame or str
        Raw APMS data in wide format (Protein column plus replicate columns
        of the form <condition>_<bait>_<rep>), or a path to a .CSV file
        containing such a table.

    Returns
    
    bait_list : list of str
        Sorted list of bait names present in the dataset.

    X_by_bait : dict
        Dictionary mapping each bait to its numeric replicate matrix
        (rows = proteins in deterministic order, columns = replicates).

    metadata : dict
        Dictionary containing:
            - "baits": the sorted bait list
            - "proteins_by_bait": dict mapping each bait to the ordered list
              of proteins used to construct the corresponding X matrix

    Notes
    
    This function performs only data loading and restructuring. It does not
    apply the EM model, tau grid search, or any statistical transformations.
    All model parameters (lambda1, lambda2, lambda3, pi, gamma) are estimated
    downstream by the hierarchical EM algorithm and are not supplied here.

    The replicate structure is preserved exactly as provided in the input.
    No averaging, collapsing, or normalization is performed at this stage.
    """

    # Case 1: user passed a DataFrame
    if isinstance(input_obj, pd.DataFrame):
        df = input_obj.copy()
    else:
        path = Path(input_obj)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {input_obj}")
        df = pd.read_csv(path)

    # Convert wide → long
    long_df = extract_bait_matrix(df)

    # Identify all baits
    bait_list = sorted(long_df["Bait"].unique())

    X_by_bait = {}
    proteins_by_bait = {}

    for bait in bait_list:
        df_b = long_df[long_df["Bait"] == bait].copy()

        # Deterministic row order
        df_b = df_b.sort_values("Protein").reset_index(drop=True)

        # Store per-bait protein list
        proteins_by_bait[bait] = df_b["Protein"].tolist()

        # Extract numeric replicate columns
        numeric_cols = [c for c in df_b.columns if c not in ("Protein", "Bait")]
        X_by_bait[bait] = df_b[numeric_cols].astype(float).to_numpy()

    metadata = {
        "baits": bait_list,
        "proteins_by_bait": proteins_by_bait,
    }

    return bait_list, X_by_bait, metadata


# %% Extract bait matrix (wide → long with replicate preservation)

def extract_bait_matrix(df):
    """
    Convert a wide-format APMS table into the long-format structure required by
    the hierarchical SAINT pipeline. The input table must contain a Protein
    column and replicate measurement columns of the form:

        <condition>_<bait>_<rep>

    where <bait> identifies the bait name and <rep> is an integer replicate
    index. Replicates are preserved exactly as provided; no collapsing or
    averaging is performed.

    This function parses the wide-format replicate columns, groups them by bait,
    and reconstructs a long-format table in which each row corresponds to a
    (Protein, Bait) pair and each replicate is placed in its own column
    (rep1, rep2, rep3, ...). Missing replicates for a bait are filled with zeros
    to maintain a consistent replicate structure.

    Parameters
    
    df : DataFrame
        Wide-format IPMS data containing a Protein column and replicate columns
        named according to the <condition>_<bait>_<rep> convention.

    Returns
    
    long_df : DataFrame
        Long-format table with columns:
            Protein : str
            Bait    : str
            rep1, rep2, rep3, ... : float
        Each bait appears once per protein, and replicate columns are ordered
        by increasing replicate index.

    Notes
    
    This function performs only structural reshaping. It does not normalize,
    filter, or transform the data. The resulting long-format table is consumed
    by `load_bait_data`, which enforces deterministic protein ordering and
    constructs the per-bait numeric matrices used by the hierarchical EM model
    and tau grid search.

    """

    protein = df["Protein"].tolist()
    value_cols = [c for c in df.columns if c != "Protein"]

    parsed = {}

    for col in value_cols:
        parts = col.split("_")
        if len(parts) < 3:
            continue

        # NAMING SCHEME:
        # <condition>_<bait>_<rep>
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

    # Build long-format dataframe
    records = []

    for bait_name, rep_dict in parsed.items():
        max_rep = max(rep_dict.keys())

        for i, p in enumerate(protein):
            row = {
                "Protein": p,
                "Bait": bait_name
            }

            for r in range(1, max_rep + 1):
                # Missing replicates → zeros
                row[f"rep{r}"] = float(rep_dict.get(r, pd.Series([0]*len(protein)))[i])

            records.append(row)

    long_df = pd.DataFrame(records)
    return long_df


