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
    Perform a two-stage grid search over the global shrinkage hyperparameter tau
    for a single bait. The function evaluates a coarse log-spaced grid followed
    by a refinement grid centered on the best coarse tau. For each tau value,
    the hierarchical EM algorithm is run to convergence and all model outputs
    and diagnostics are recorded.

    Parameters
    
    X : array-like, shape (n_proteins, n_replicates)
        The data matrix for the current bait. Each row corresponds to a protein
        and each column corresponds to a replicate measurement.

    hyperparams : dict
        Hyperparameters for the hierarchical EM model. Any missing entries are
        filled internally by the EM wrapper. Only the "tau" entry is modified
        during the grid search; all other hyperparameters remain unchanged.

    biological_bait : str
        Identifier for the current bait. Used only for labeling and diagnostics.

    Returns
    
    dict
        A dictionary containing the following entries:

        best_tau : float
            The tau value that produced the highest final log-likelihood across
            both the coarse and refinement grids.

        best_result : dict
            The full EM output corresponding to best_tau. This includes:
                - lambda1, lambda2, lambda3 : component means
                - pi : mixture proportions
                - gamma : posterior component probabilities
                - loglik_history : list of log-likelihood values per EM iteration
                - convergence_info : dict describing EM convergence diagnostics
                - iteration_count : number of EM iterations performed

        tau_grid_results : dict
            A mapping from each tau value (float) to the full EM result for that
            tau. Each entry contains the same fields as best_result.

        convergence_info : dict
            A mapping from each tau value to the EM convergence diagnostics
            returned by the EM wrapper. This allows inspection of EM behavior
            across the entire grid.

        iteration_counts : dict
            A mapping from each tau value to the number of EM iterations
            performed for that tau.

        tau_grid : list of float
            The complete set of tau values evaluated across both the coarse and
            refinement grids, sorted in ascending order.

    Notes
    
    - The coarse grid is log-spaced between 10^COARSE_MIN and 10^COARSE_MAX.
    - The refinement grid is log-spaced around log10(best_tau) with a half-range
      of REFINE_HALF_RANGE.
    - Only the tau entry in the hyperparameter dictionary is modified during the
      search. All other hyperparameters remain unchanged.
    - The EM wrapper must return loglik_history, convergence_info, and
      iteration_count for diagnostics to be recorded.
    """

    # %% Stage 1: coarse grid
    coarse_grid = np.logspace(COARSE_MIN, COARSE_MAX, num=COARSE_POINTS)

    tau_grid_results = {}
    convergence_info = {}
    iteration_counts = {}

    best_tau = None
    best_loglik = -np.inf
    best_result = None

    for tau in coarse_grid:
        hyper = dict(hyperparams)
        hyper["tau"] = float(tau)

        # FIX: use hyper, not hyperparams
        result = run_em_hierarchical(
            X=X,
            hyperparams=hyper,
            biological_bait=biological_bait,
            max_iter=hyper.get("max_iter", 100),
            tol_loglik=hyper.get("tol_loglik", 1e-6),
            tol_params=hyper.get("tol_params", 1e-6),
            seed=hyper.get("seed", None),
            )
        tau_grid_results[tau] = result

        # Extract diagnostics
        convergence_info[tau] = result.get("convergence_info", {})
        iteration_counts[tau] = result.get("iteration_count", None)

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

        # Extract diagnostics
        convergence_info[tau] = result.get("convergence_info", {})
        iteration_counts[tau] = result.get("iteration_count", None)

        final_loglik = result["loglik_history"][-1]
        if final_loglik > best_loglik:
            best_loglik = final_loglik
            best_tau = tau
            best_result = result

    # %% Final output
    return {
        "best_tau": best_tau,
        "best_result": best_result,
        "tau_grid_results": tau_grid_results,
        "convergence_info": convergence_info,
        "iteration_counts": iteration_counts,
        "tau_grid": sorted(tau_grid_results.keys()),
    }