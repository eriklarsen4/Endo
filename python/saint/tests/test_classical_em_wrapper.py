# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:57:43 2026

@author: Erik
"""

# %% Imports
import numpy as np
from saint.model.classical_em_wrapper import run_em_classical


# %% Tests

def test_em_runs_and_returns_expected_keys():
    # X_sum must be a 1D vector because run_em_classical collapses counts
    X_sum = np.array([5.0, 3.0])

    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([2.0, 2.0]),
        "pi_init": np.array([0.7, 0.3]),
        "alpha": np.array([2.0, 2.0]),   # unused but harmless
    }

    results = run_em_classical(
        X_sum,
        hyperparams,
        "BAIT",
        max_iter=5,
        seed=1,
    )

    expected_keys = {
        "loglik_history",
        "lambda1_history",
        "lambda2_history",
        "pi_history",
        "gamma_history",
        "lambda1",
        "lambda2",
        "pi",
        "gamma",
    }

    assert expected_keys.issubset(results.keys())


def test_em_loglik_increases_or_stabilizes():
    X_sum = np.array([4.0, 4.0])

    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([3.0, 3.0]),
        "pi_init": np.array([0.5, 0.5]),
        "alpha": np.array([2.0, 2.0]),
    }

    results = run_em_classical(
        X_sum,
        hyperparams,
        "BAIT",
        max_iter=10,
        seed=0,
    )

    loglik = results["loglik_history"]

    # EM log-likelihood must be non-decreasing
    assert all(loglik[i] <= loglik[i+1] + 1e-8 for i in range(len(loglik) - 1))
