# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:39:06 2026

@author: Erik
"""

# %% Imports
import pandas as pd
import numpy as np

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline


# %% Tests

def test_tau_grid_integration_structure():
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

    # Run the full hierarchical pipeline (tau-grid included)
    results = run_hierarchical_pipeline(
        input_data,
        bait_names,
        metadata,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    raw = results["raw_outputs"]

    # Tau info must exist for each bait
    assert "tau_info" in raw
    assert isinstance(raw["tau_info"], dict)
    assert "B1" in raw["tau_info"]

    tau_grid_result = raw["tau_info"]["B1"]

    # Expected tau-grid keys
    expected_keys = {
        "taus",
        "logliks",
        "em_results",
        "convergence_info",
        "iteration_counts",
        "tolerance",
    }
    assert expected_keys.issubset(tau_grid_result.keys())

    # Basic structural checks
    taus = tau_grid_result["taus"]
    logliks = tau_grid_result["logliks"]
    em_results = tau_grid_result["em_results"]

    assert isinstance(taus, (list, np.ndarray))
    assert len(taus) == len(logliks)
    assert set(em_results.keys()) == set(taus)

    # Check EM result structure for one tau
    tau0 = taus[0]
    em0 = em_results[tau0]

    for key in [
        "lambda1", "lambda2", "lambda3",
        "pi", "gamma",
        "loglik_history",
        "lambda1_history", "lambda2_history", "lambda3_history",
        "pi_history", "gamma_history",
        "alpha_history", "a_history", "b_history",
        "convergence_info", "iteration_count",
    ]:
        assert key in em0