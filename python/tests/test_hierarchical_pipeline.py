# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:30:57 2026

@author: Erik
"""

#%% Import
import pandas as pd

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline

#%% tests
def test_hierarchical_pipeline_runs_and_returns_expected_structure():
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

    results = run_hierarchical_pipeline(
        input_data=input_data,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    # Top-level keys
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata structure (dict, with new fields)
    meta = results["metadata"]
    assert isinstance(meta, dict)
    for key in ["baits", "proteins_by_bait", "control_baits", "treatment_baits"]:
        assert key in meta

    # Raw outputs structure
    raw = results["raw_outputs"]
    assert isinstance(raw, dict)
    assert "em_results" in raw
    assert "tau_info" in raw

    # EM results must contain all baits listed in metadata
    assert set(raw["em_results"].keys()) == set(meta["baits"])

    # Results table
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Core columns: protein ID, bait ID, and hierarchical gamma3
    for col in ["Protein", "bait", "gamma3"]:
        assert col in df.columns