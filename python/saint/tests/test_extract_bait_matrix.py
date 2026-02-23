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
    # Wide-format synthetic input
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

    # Expected columns
    expected_cols = {"Protein", "Bait", "rep1", "rep2"}
    assert expected_cols.issubset(long_df.columns)

    # Should have 2 proteins × 2 baits = 4 rows
    assert len(long_df) == 4

    # Check one bait block explicitly
    b1_rows = long_df[long_df["Bait"] == "B1"].sort_values("Protein")
    assert np.allclose(b1_rows["rep1"].to_numpy(), [5.0, 1.0])
    assert np.allclose(b1_rows["rep2"].to_numpy(), [3.0, 0.0])

    # Check control bait block
    ctrl_rows = long_df[long_df["Bait"] == "CTRL"].sort_values("Protein")
    assert np.allclose(ctrl_rows["rep1"].to_numpy(), [2.0, 0.0])
    assert np.allclose(ctrl_rows["rep2"].to_numpy(), [1.0, 1.0])


def test_extract_bait_matrix_missing_replicates_filled_with_zero():
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condA_B1_1": [5, 1],
        # B1 has no rep2 column → should be filled with zeros
        "condA_CTRL_1": [2, 0],
        "condA_CTRL_2": [1, 1],
    })

    long_df = extract_bait_matrix(df)

    # B1 should have rep1 and rep2, with rep2 filled with zeros
    b1_rows = long_df[long_df["Bait"] == "B1"].sort_values("Protein")

    assert np.allclose(b1_rows["rep1"].to_numpy(), [5.0, 1.0])
    assert np.allclose(b1_rows["rep2"].to_numpy(), [0.0, 0.0])