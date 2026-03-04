# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:25:10 2026

@author: Erik
"""
# %% Import
import numpy as np
import pandas as pd

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline

# %% test
def test_pipeline_core_structure(minimal_input_data):
    """Pipeline should run end-to-end and return the correct top-level keys."""

    out = run_hierarchical_pipeline(
        input_data=minimal_input_data,
        make_plots=False,
        verbose=False,
    )

    assert "raw_outputs" in out
    assert "metadata" in out
    assert "results_df" in out

    # Check EM results structure
    em_results = out["raw_outputs"]["em_results"]
    assert isinstance(em_results, dict)
    assert len(em_results) > 0

    # Check results_df schema
    df = out["results_df"]
    expected_cols = [
        "Protein", "bait",
        "lambda1", "lambda2", "lambda3",
        "tau",
        "pi1", "pi2", "pi3",
        "gamma1", "gamma2", "gamma3",
    ]
    assert list(df.columns) == expected_cols
    assert len(df) > 0