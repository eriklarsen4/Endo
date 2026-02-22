# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 18:22:55 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd

from saint.model.hierarchical_em_wrapper import run_em_hierarchical
from saint.io.data_input import extract_bait_matrix


# %% Hierarchical SAINT pipeline

def run_hierarchical_pipeline(
    input_data,
    bait_names,
    metadata,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False,
    make_plots=False,
    plot_dir="plots_hierarchical"
):
    """
    Run the hierarchical SAINT pipeline. This function preserves replicate level counts,
    runs hierarchical EM for each bait, merges final parameter values into the original
    input dataframe, and returns a unified results object containing raw EM histories,
    metadata, and a sorted results dataframe.

    input_data is the original wide format dataframe.
    bait_names is the list of experimental baits.
    metadata contains biological bait names, accession numbers, and molecular weights.
    max_iter, tol_loglik, tol_params, seed, verbose, make_plots, and plot_dir control EM behavior.

    The output is a dictionary with keys raw_outputs, metadata, and results_df.
    """

    # %% Extract long format data (preserves replicate columns)
    long_df = extract_bait_matrix(input_data)

    # %% Identify controls
    all_baits = sorted(long_df["Bait"].unique())
    controls = [b for b in all_baits if b not in bait_names]

    controls_used = {}
    for ctrl in controls:
        ctrl_rows = long_df[long_df["Bait"] == ctrl]
        rep_cols_ctrl = [c for c in ctrl_rows.columns if c.startswith("rep")]
        controls_used[ctrl] = rep_cols_ctrl

    # %% Prepare output containers
    raw_outputs = {}
    results_rows = []

    # %% Run EM per bait
    for bait in bait_names:

        df_bait = long_df[long_df["Bait"] == bait].copy()

        rep_cols = [c for c in df_bait.columns if c.startswith("rep")]
        X = df_bait[rep_cols].to_numpy()

        biological_bait = metadata["biological_bait_names"][bait]

        X_sum = X.sum(axis=1).astype(float)
        mean_level = max(X_sum.mean(), 1.0)

        hyperparams = {
            "lambda1_init": np.full(X.shape[0], 0.5 * mean_level),
            "lambda2_init": np.full(X.shape[0], 1.0 * mean_level),
            "lambda3_init": np.full(X.shape[0], 2.0 * mean_level),
            "pi_init": np.array([0.6, 0.3, 0.1], dtype=float)
        }

        results_em = run_em_hierarchical(
            X,
            hyperparams,
            biological_bait,
            max_iter=max_iter,
            tol_loglik=tol_loglik,
            tol_params=tol_params,
            seed=seed,
            verbose=verbose
        )

        raw_outputs[bait] = {
            "loglik": results_em["loglik_history"],
            "lambda1": results_em["lambda1_history"],
            "lambda2": results_em["lambda2_history"],
            "lambda3": results_em["lambda3_history"],
            "tau": results_em["tau_history"],
            "pi": results_em["pi_history"],
            "gamma": results_em["gamma_history"],
            "alpha": results_em["alpha_history"],
            "a": results_em["a_history"],
            "b": results_em["b_history"]
        }

        lambda1 = results_em["lambda1"]
        lambda2 = results_em["lambda2"]
        lambda3 = results_em["lambda3"]
        tau = results_em["tau"]
        pi = results_em["pi"]
        gamma = results_em["gamma"]

        df_bait_reset = df_bait.reset_index(drop=True)

        for i, row in df_bait_reset.iterrows():
            results_rows.append({
                "Protein": row["Protein"],
                "Bait": bait,
                **{col: row[col] for col in rep_cols},
                "lambda1": lambda1[i],
                "lambda2": lambda2[i],
                "lambda3": lambda3[i],
                "tau": tau[i],
                "pi1": pi[0],
                "pi2": pi[1],
                "pi3": pi[2],
                "gamma1": gamma[i, 0],
                "gamma2": gamma[i, 1],
                "gamma3": gamma[i, 2]
            })

    # %% Build results_df
    results_df = pd.DataFrame(results_rows)

    # %% Sort by Protein then gamma3 descending
    results_df = results_df.sort_values(
        by=["Protein", "gamma3"],
        ascending=[True, False],
        ignore_index=True
    )

    # %% Build unified metadata
    metadata_out = {
        "bait_names": bait_names,
        "biological_bait_names": metadata["biological_bait_names"],
        "AN": metadata["AN"],
        "MW": metadata["MW"],
        "controls_used": controls_used
    }

    # %% Final unified output
    return {
        "raw_outputs": raw_outputs,
        "metadata": metadata_out,
        "results_df": results_df
    }

