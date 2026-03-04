# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:37:51 2026

@author: Erik
"""

# %% Imports
import numpy as np
from saint.model.hierarchical_tau_grid import run_tau_grid


# %% Tests

def test_tau_grid_runs_and_returns_expected_keys():
    # Minimal synthetic data
    X = np.array([[5, 3], [1, 0]])

    # Minimal hyperparameters (only those used by EM)
    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([2.0, 2.0]),
        "lambda3_init": np.array([4.0, 4.0]),
        "pi_init": np.array([0.5, 0.3, 0.2]),
    }

    # A tiny tau grid
    tau_values = [0.1, 0.5]

    results = run_tau_grid(
        X=X,
        hyperparams=hyperparams,
        biological_bait="BAIT",
        max_iter=5,
        seed=1,
    )

    expected_keys = {
        "taus",
        "logliks",
        "em_results",
        "convergence_info",
        "iteration_counts",
        "tolerance",
    }

    assert expected_keys.issubset(results.keys())

    # Structural checks
    assert list(results["taus"]) == tau_values
    assert len(results["logliks"]) == len(tau_values)
    assert set(results["em_results"].keys()) == set(tau_values)

    # EM result structure sanity check
    em0 = results["em_results"][tau_values[0]]
    assert "lambda1" in em0
    assert "lambda2" in em0
    assert "lambda3" in em0
    assert "pi" in em0
    assert "gamma" in em0
    assert "loglik_history" in em0