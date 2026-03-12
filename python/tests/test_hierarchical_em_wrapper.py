# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:34:40 2026

@author: Erik
"""

# %% Imports
import numpy as np
from saint.model.hierarchical_em_wrapper import run_em_hierarchical


# %% Tests

def test_em_hierarchical_runs_and_returns_expected_keys():
    X = np.array([[5, 3], [1, 0]])

    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([2.0, 2.0]),
        "lambda3_init": np.array([4.0, 4.0]),
        "pi_init": np.array([0.5, 0.3, 0.2]),
    }

    results = run_em_hierarchical(
        X=X,
        hyperparams=hyperparams,
        biological_bait="BAIT",
        max_iter=5,
        seed=1,
    )

    expected_keys = {
        "loglik_history",
        "lambda1_history",
        "lambda2_history",
        "lambda3_history",
        "pi_history",
        "gamma_history",
        "alpha_history",
        "a_history",
        "b_history",
        "lambda1",
        "lambda2",
        "lambda3",
        "tau",
        "pi",
        "gamma",
    }

    assert expected_keys.issubset(results.keys())

