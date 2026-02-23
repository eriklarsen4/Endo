# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:08:51 2026

@author: Erik
"""

# %% Imports
import numpy as np
import pandas as pd

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline


# %% Tests

def test_hierarchical_pipeline_runs_and_returns_expected_structure():
    # Minimal synthetic wide-format input
    input_data = pd.DataFrame({
        "Protein": ["P1", "P2", "P1", "P2"],
        "Bait":    ["B1", "B1", "CTRL", "CTRL"],
        "rep1":    [5, 1, 2, 0],
        "rep2":    [3, 0, 1, 1],
    })

    bait_names = ["B1"]

    metadata = {
        "biological_bait_names": {"B1": "BAIT1"},
        "AN": {"B1": "AN1"},
        "MW": {"B1": 50.0},
    }

    results = run_hierarchical_pipeline(
        input_data,
        bait_names,
        metadata,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    # Top-level keys
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata structure
    meta = results["metadata"]
    assert meta["bait_names"] == bait_names
    assert "biological_bait_names" in meta
    assert "AN" in meta
    assert "MW" in meta
    assert "controls_used" in meta

    # Raw outputs
    raw = results["raw_outputs"]
    assert "B1" in raw

    # EM histories must exist
    expected_raw_keys = {
        "loglik", "lambda1", "lambda2", "lambda3",
        "tau", "pi", "gamma",
        "alpha", "a", "b",
    }
    assert expected_raw_keys.issubset(raw["B1"].keys())

    # Results dataframe
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)

    expected_cols = {
        "Protein", "Bait", "rep1", "rep2",
        "lambda1", "lambda2", "lambda3",
        "tau",
        "pi1", "pi2", "pi3",
        "gamma1", "gamma2", "gamma3",
    }
    assert expected_cols.issubset(df.columns)

    # Sorting check: Protein ascending, gamma3 descending
    df_sorted = df.sort_values(
        by=["Protein", "gamma3"],
        ascending=[True, False],
        ignore_index=True,
    )
    assert df.equals(df_sorted)