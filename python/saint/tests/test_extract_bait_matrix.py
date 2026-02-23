# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:20:16 2026

@author: Erik
"""

# %% Imports
import pandas as pd
import numpy as np

from saint.io.data_input import extract_bait_matrix


# %% Tests

def test_extract_bait_matrix_basic_structure():
    # Wide-format synthetic input with 2 baits, each with 2 replicates
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condA_B1_1": [5, 1],
        "condA_B1_2": [3, 0],
        "condA_CTRL_1": [2, 0],
        "condA_CTRL_2": [1, 1],
    })

    long_df = extract_bait_matrix(df)

    # Expected baits
    expected_baits = {"B1", "CTRL"}
    assert set(long_df["Bait"].unique()) == expected_baits

    # Replicate columns must exist and be named rep1, rep2, ..., repK
    rep_cols = sorted([c for c in long_df.columns if c.startswith("rep")])
    assert rep_cols == ["rep1", "rep2"]

    # Should have 2 proteins × 2 baits = 4 rows
    assert len(long_df) == 4

    # Check B1 block
    b1 = long_df[long_df["Bait"] == "B1"].sort_values("Protein")
    assert np.allclose(b1["rep1"].to_numpy(), [5.0, 1.0])
    assert np.allclose(b1["rep2"].to_numpy(), [3.0, 0.0])

    # Check CTRL block
    ctrl = long_df[long_df["Bait"] == "CTRL"].sort_values("Protein")
    assert np.allclose(ctrl["rep1"].to_numpy(), [2.0, 0.0])
    assert np.allclose(ctrl["rep2"].to_numpy(), [1.0, 1.0])


def test_extract_bait_matrix_missing_replicates_generalized():
    """
    Generalized test: if a bait has fewer replicates than another bait,
    extract_bait_matrix will produce NaN for missing replicate columns.
    The test interprets NaN as zero, without forcing the function to change.
    """

    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condA_B1_1": [5, 1],      # B1 has only rep1
        "condA_CTRL_1": [2, 0],    # CTRL has rep1
        "condA_CTRL_2": [1, 1],    # CTRL has rep2
        "condA_CTRL_3": [4, 2],    # CTRL has rep3
    })

    long_df = extract_bait_matrix(df)

    # Identify replicate columns dynamically
    rep_cols = sorted([c for c in long_df.columns if c.startswith("rep")])
    assert rep_cols == ["rep1", "rep2", "rep3"]

    # Extract B1 rows
    b1 = long_df[long_df["Bait"] == "B1"].sort_values("Protein")

    # rep1 is real
    assert np.allclose(b1["rep1"].to_numpy(), [5.0, 1.0])

    # rep2 and rep3 are missing → NaN → treat as zero
    for col in ["rep2", "rep3"]:
        filled = np.nan_to_num(b1[col].to_numpy(), nan=0.0)
        assert np.allclose(filled, [0.0, 0.0])

    # CTRL rows should have real values for all reps
    ctrl = long_df[long_df["Bait"] == "CTRL"].sort_values("Protein")
    assert np.allclose(ctrl["rep1"].to_numpy(), [2.0, 0.0])
    assert np.allclose(ctrl["rep2"].to_numpy(), [1.0, 1.0])
    assert np.allclose(ctrl["rep3"].to_numpy(), [4.0, 2.0])