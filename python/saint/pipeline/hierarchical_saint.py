# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 18:22:55 2026

@author: Erik
"""

# %% Imports
import numpy as np
import pandas as pd
from saint.io.data_input import load_bait_data
from saint.model.hierarchical_tau_grid import run_tau_grid


# %% Hierarchical pipeline
def run_hierarchical_pipeline(
    input_data,
    hyperparams=None,
    metadata=None,
    bait_names=None,
    make_plots=False,
    plot_dir=None
):
    """
    Run the hierarchical SAINT pipeline for all baits using a per bait tau grid
    search. This function loads the data, applies the hierarchical EM model with
    a two stage tau grid for each bait, collects the best fitting results, and
    assembles a unified prey level results table. The hyperparameter dictionary
    is passed to the EM wrapper and defaults are filled internally if not
    supplied. Only the tau entry is modified during the grid search.
    
    Parameters
    
    input_data : dict or DataFrame
        Raw APMS data in the format expected by the data loader.
    hyperparams : dict
        Optional hyperparameters for the hierarchical EM model. Users may supply
        overrides but defaults are constructed internally for any missing entries.
        The tau entry is overridden during the grid search.
        
    Returns
    
    dict
        Unified results object containing raw outputs from the hierarchical EM
        model with tau tuning, metadata describing the experiment, and a combined
        prey level results DataFrame.
        
   Model parameters
   
   The model parameters lambda1, lambda2, lambda3, pi, and gamma are estimated
   by the EM algorithm and are not supplied by the user.
   
   Hyperparameters
   
   The hyperparameters alpha, a_k, b_k, and tau may be supplied in the
   hyperparameter dictionary but are not required. Any missing entries are filled
   with internal defaults. The tau entry is overridden by the grid search. All
   other hyperparameters remain unchanged unless explicitly provided.
   
   Tau grid tuning parameters
   
   The coarse and refinement grid ranges and resolutions define the search over
   tau. These are tuning parameters for the grid search procedure and are not
   model hyperparameters.
    """

    # %% Load data
    bait_list, X_by_bait, metadata = load_bait_data(input_data)

    # %% Storage
    all_results = {}
    all_tau_info = {}

    # %% Parameter initialization
    if hyperparams is None:
        hyperparams = {}

    # Dirichlet prior (alpha = 2 for each component)
    hyperparams.setdefault("alpha1", 2.0)
    hyperparams.setdefault("alpha2", 2.0)
    hyperparams.setdefault("alpha3", 2.0)
    
    # Mixing proportions initialization (uniform)
    hyperparams.setdefault("pi_init", np.ones(3, dtype=float) / 3.0)
    
    # Tau grid (your actual coarse grid: logspace(-1, 1, 7))
    hyperparams.setdefault("tau_grid", np.logspace(-1.0, 1.0, num=7))
    
    # Tau default = geometric center of the grid = 10^0 = 1.0
    hyperparams.setdefault("tau", 1.0)

    # %% Per bait tau grid
    for bait in bait_list:
        X = X_by_bait[bait]

        tau_output = run_tau_grid(X, hyperparams, bait)

        all_results[bait] = tau_output["best_result"]
        all_tau_info[bait] = {
            "best_tau": tau_output["best_tau"],
            "tau_grid_results": tau_output["tau_grid_results"]
        }

    # %% Build results_df
    results_df = _assemble_results_df(all_results, metadata)

    # %% Final unified output
    return {
        "raw_outputs": {
            "em_results": all_results,
            "tau_info": all_tau_info
        },
        "metadata": metadata,
        "results_df": results_df
    }


# %% Helper to assemble results_df
def _assemble_results_df(all_results, metadata):
    """
    Assemble a single results dataframe from per-bait EM outputs.
    Expects each result dict to contain lambda1, lambda2, lambda3, tau, pi, gamma.
    Uses metadata["proteins"] for the prey names.
    """
    proteins = metadata["proteins"]
    rows = []

    for bait, result in all_results.items():
        proteins = metadata["proteins_by_bait"][bait]
        lambda1 = result["lambda1"]
        lambda2 = result["lambda2"]
        lambda3 = result["lambda3"]
        tau = result["tau"]
        pi = result["pi"]
        gamma = result["gamma"]

        df = pd.DataFrame({
            "Protein": proteins,
            "bait": bait,
            "lambda1": lambda1,
            "lambda2": lambda2,
            "lambda3": lambda3,
            "tau": tau,
            "pi1": pi[0],
            "pi2": pi[1],
            "pi3": pi[2],
            "gamma1": gamma[:, 0],
            "gamma2": gamma[:, 1],
            "gamma3": gamma[:, 2],
        })

        rows.append(df)

    results_df = pd.concat(rows, ignore_index=True)
    results_df = results_df.sort_values(["Protein", "gamma3"], ascending=[True, False])
    return results_df

