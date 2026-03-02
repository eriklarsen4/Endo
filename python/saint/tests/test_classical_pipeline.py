# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:07:19 2026

@author: Erik
"""

# %% Imports
import numpy as np
import pandas as pd

from saint.pipeline.classical_saint import run_classical_pipeline


# %% Tests

def test_classical_pipeline_runs_and_returns_expected_structure():
    # Minimal synthetic wide-format input
    input_data = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condA_B1_1": [5, 1],
        "condA_B1_2": [3, 0],
        "condA_CTRL_1": [2, 0],
        "condA_CTRL_2": [1, 1],
    })

    bait_names = ["B1"]

    metadata = {
        "biological_bait_names": {"B1": "BAIT1"},
        "AN": {"B1": "AN1"},
        "MW": {"B1": 50.0},
    }

    results = run_classical_pipeline(
        input_data,
        bait_names,
        metadata,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    # Top-level keys
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata object (typed)
    meta = results["metadata"]
    assert hasattr(meta, "user_provided_fields")
    assert hasattr(meta, "inferred_fields")
    assert hasattr(meta, "pipeline_derived_fields")

    # Pipeline-derived fields
    pdf = meta.pipeline_derived_fields
    assert pdf.bait_names == bait_names
    assert isinstance(pdf.hyperparameters, dict)
    assert pdf.tau_grid == []  # classical pipeline has no tau grid
    assert isinstance(pdf.convergence, dict)
    assert isinstance(pdf.iteration_counts, dict)

    # Raw outputs structure
    raw = results["raw_outputs"]
    assert set(raw.keys()) == {"em_results", "tau_info"}

    # EM results must contain the bait
    assert "B1" in raw["em_results"]

    em = raw["em_results"]["B1"]

    # EM histories and final values must exist
    for key in [
        "loglik_history",
        "lambda1_history",
        "lambda2_history",
        "pi_history",
        "gamma_history",
        "lambda1",
        "lambda2",
        "pi",
        "gamma",
        "convergence_info",
        "iteration_count",
    ]:
        assert key in em

    # Results dataframe
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)

    expected_cols = {
        "Protein",
        "bait",
        "lambda1",
        "lambda2",
        "lambda3",
        "tau",
        "pi1",
        "pi2",
        "pi3",
        "gamma1",
        "gamma2",
        "gamma3",
    }
    assert expected_cols.issubset(df.columns)

    # Sorting check: Protein ascending, gamma3 descending
    df_sorted = df.sort_values(
        by=["Protein", "gamma3"],
        ascending=[True, False],
        ignore_index=True,
    )
    assert df.equals(df_sorted)