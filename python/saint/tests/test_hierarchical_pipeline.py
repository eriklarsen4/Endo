# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:15:20 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd
from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline


# %% Tests

def test_hierarchical_pipeline_runs_and_returns_expected_structure():
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "Bait": ["B1", "B1"],
        "Count_1": [5, 1],
        "Count_2": [3, 0]
    })

    metadata = {"bait_to_biological_name": {"B1": "B1"}}

    results = run_hierarchical_pipeline(
        input_data=df,
        bait_names=["B1"],
        metadata=metadata,
        max_iter=5,
        seed=0,
        make_plots=False
    )

    assert "B1" in results
    entry = results["B1"]

    expected_keys = {"output_df", "merged_df", "raw_outputs", "metadata"}
    assert expected_keys.issubset(entry.keys())

    out = entry["output_df"]
    assert set(out.columns) == {
        "Protein",
        "lambda1",
        "lambda2",
        "lambda3",
        "tau",
        "pi1",
        "pi2",
        "pi3",
        "gamma1",
        "gamma2",
        "gamma3"
    }

