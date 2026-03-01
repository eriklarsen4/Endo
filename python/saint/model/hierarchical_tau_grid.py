# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 09:02:22 2026

@author: Erik
"""

# %% Imports
import numpy as np
from saint.model.hierarchical_em_wrapper import run_em_hierarchical


# %% Tau grid parameters
COARSE_MIN = -1.0
COARSE_MAX = 1.0
COARSE_POINTS = 7

REFINE_HALF_RANGE = 0.5
REFINE_POINTS = 7


# %% Tau grid search
def run_tau_grid(X, hyperparams, biological_bait):
    """
    Perform a two stage grid search for the global shrinkage parameter tau.
    The function evaluates a coarse log spaced grid followed by a refinement
    grid centered on the best coarse tau. Full EM results are stored for each
    tau value. Only the tau entry in the hyperparameter dictionary is modified.
    All other hyperparameters are filled by the EM wrapper if not provided.

    Parameters
    ----------
    X : array
        Data matrix for the current bait.
    hyperparams : dict
        Hyperparameters for the hierarchical EM model. Users may supply values
        but defaults are filled inside the EM wrapper. Only the tau entry is
        overridden during the grid search.
    biological_bait : str
        Identifier for the current bait.

    Returns
    -------
    dict
        Dictionary containing best_tau, best_result, and tau_grid_results.
    """

    # %% Stage 1: coarse grid
    coarse_grid = np.logspace(COARSE_MIN, COARSE_MAX, num=COARSE_POINTS)
    tau_grid_results = {}
    best_tau = None
    best_loglik = -np.inf
    best_result = None

    for tau in coarse_grid:
        hyper = dict(hyperparams)
        hyper["tau"] = float(tau)

        result = run_em_hierarchical(X, hyper, biological_bait)
        tau_grid_results[tau] = result

        final_loglik = result["loglik_history"][-1]
        if final_loglik > best_loglik:
            best_loglik = final_loglik
            best_tau = tau
            best_result = result

    # %% Stage 2: refinement grid
    center = np.log10(best_tau)
    refine_grid = np.logspace(center - REFINE_HALF_RANGE,
                              center + REFINE_HALF_RANGE,
                              num=REFINE_POINTS)

    for tau in refine_grid:
        hyper = dict(hyperparams)
        hyper["tau"] = float(tau)

        result = run_em_hierarchical(X, hyper, biological_bait)
        tau_grid_results[tau] = result

        final_loglik = result["loglik_history"][-1]
        if final_loglik > best_loglik:
            best_loglik = final_loglik
            best_tau = tau
            best_result = result

    # %% Final output
    return {
        "best_tau": best_tau,
        "best_result": best_result,
        "tau_grid_results": tau_grid_results
    }