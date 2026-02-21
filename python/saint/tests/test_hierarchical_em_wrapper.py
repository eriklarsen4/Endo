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
        "tau_init": np.array([1.0, 1.0]),
        "pi_init": np.array([0.5, 0.3, 0.2]),
        "alpha": np.array([2.0, 2.0, 2.0])
    }

    results = run_em_hierarchical(
        X=X,
        hyperparams=hyperparams,
        biological_bait="BAIT",
        max_iter=5,
        seed=1
    )

    expected_keys = {
        "lambda1",
        "lambda2",
        "lambda3",
        "tau",
        "pi",
        "gamma",
        "loglik_history",
        "lambda1_history",
        "lambda2_history",
        "lambda3_history",
        "tau_history",
        "pi_history",
        "gamma_history",
        "biological_bait"
    }

    assert expected_keys.issubset(results.keys())


def test_em_hierarchical_loglik_increases_or_stabilizes():
    X = np.array([[4, 4], [2, 2]])
    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([3.0, 3.0]),
        "lambda3_init": np.array([6.0, 6.0]),
        "tau_init": np.array([1.0, 1.0]),
        "pi_init": np.array([0.4, 0.4, 0.2]),
        "alpha": np.array([2.0, 2.0, 2.0])
    }

    results = run_em_hierarchical(
        X=X,
        hyperparams=hyperparams,
        biological_bait="BAIT",
        max_iter=10,
        seed=0
    )

    loglik = results["loglik_history"]
    assert all(loglik[i] <= loglik[i+1] + 1e-8 for i in range(len(loglik)-1))

