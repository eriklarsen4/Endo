# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:06:37 2026

@author: Erik
"""

# %% Imports

import pandas as pd
from saint.io.data_input import extract_bait_matrix


# %% Test reshape preserves replicate columns

def test_extract_bait_matrix_rep_preservation():
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condX_BaitA_1": [5, 1],
        "condX_BaitA_2": [7, 0],
        "condX_BaitB_1": [3, 0],
        "condX_BaitB_2": [4, 1]
    })

    long_df = extract_bait_matrix(df)

    assert set(long_df.columns) == {
        "Protein", "Bait", "rep1", "rep2"
    }

    df_a = long_df[long_df["Bait"] == "BaitA"]
    assert df_a["rep1"].tolist() == [5.0, 1.0]
    assert df_a["rep2"].tolist() == [7.0, 0.0]

    df_b = long_df[long_df["Bait"] == "BaitB"]
    assert df_b["rep1"].tolist() == [3.0, 0.0]
    assert df_b["rep2"].tolist() == [4.0, 1.0]

