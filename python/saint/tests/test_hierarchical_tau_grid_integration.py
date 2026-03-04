# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:39:06 2026

@author: Erik
"""

# %% Import
import pandas as pd

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline

# %% Tests
def test_tau_grid_integration_structure():
    # Minimal synthetic wide-format input
    input_data = pd.DataFrame(
        {
            "Protein": ["P1", "P2"],
            "condA_B1_1": [5, 1],
            "condA_B1_2": [3, 0],
            "condA_CTRL_1": [2, 0],
            "condA_CTRL_2": [1, 1],
        }
    )

    # Run the full hierarchical pipeline (tau-grid included)
    results = run_hierarchical_pipeline(
        input_data=input_data,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    # Top-level structure
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata structure
    meta = results["metadata"]
    assert isinstance(meta, dict)
    for key in ["baits", "proteins_by_bait", "control_baits", "treatment_baits"]:
        assert key in meta

    # Tau-grid integration
    raw = results["raw_outputs"]
    assert "tau_info" in raw
    tau_info = raw["tau_info"]

    # tau_info may be empty or populated depending on pipeline settings,
    # but it must be a dict.
    assert isinstance(tau_info, dict)

    # EM results must exist for each bait
    assert "em_results" in raw
    assert set(raw["em_results"].keys()) == set(meta["baits"])

    # Results table
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Protein" in df.columns
    assert "bait" in df.columns
    assert "gamma3" in df.columns