# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 18:22:55 2026

@author: Erik
"""

# %% Imports
from saint.io.data_input import load_bait_data
from saint.model.hierarchical_tau_grid import run_tau_grid


# %% Hierarchical pipeline
def run_hierarchical_pipeline(input_data, hyperparams):
    """
    Run the hierarchical SAINT pipeline for all baits using a per bait tau grid search.
    This function loads the data, applies the hierarchical EM model with a two-stage
    tau grid for each bait, collects the best fitting results, and assembles a unified
    prey-level results table. The hyperparameter dictionary is passed to the EM wrapper
    and defaults are filled internally if not supplied. Only the tau entry is modified
    during the grid search.

    Parameters
    ----------
    input_data : dict or DataFrame
        Raw APMS data in the format expected by the data loader.
    hyperparams : dict
        Hyperparameters for the hierarchical EM model. Users may supply values but
        defaults are filled inside the EM wrapper. The tau entry is overridden during
        the grid search.

    Returns
    -------
    dict
        Unified results object containing raw outputs from the hierarchical EM model
        with tau tuning, metadata describing the experiment, and a combined prey level
        results DataFrame.
    """

    # %% Load data
    bait_list, X_by_bait, metadata = load_bait_data(input_data)

    # %% Storage
    all_results = {}
    all_tau_info = {}

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
    Convert per bait EM results into a unified prey level results DataFrame.
    """
    rows = []

    for bait, result in all_results.items():
        df = result["results_df"].copy()
        df["bait"] = bait
        rows.append(df)

    return (
        rows[0].__class__.concat(rows, ignore_index=True)
        if rows else rows
    )

