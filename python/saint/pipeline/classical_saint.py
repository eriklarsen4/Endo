# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 23:41:42 2026

@author: Erik
"""

# %% Imports

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from saint.io.data_input import extract_bait_matrix
from saint.model.classical_em_wrapper import run_em
from saint.diagnostics.diagnostics_classical import make_classical_plots


# %% Helpers

def reshape_counts_df(df, bait, biological_bait):
    """
    Extract the count matrix X for a specific bait.

    X is a numeric matrix of prey (protein) by experiment.
    df_bait is the subset of the input DataFrame corresponding to the bait.
    """
    df_bait = df[df["Bait"] == bait].copy()
    count_cols = [c for c in df_bait.columns if c.startswith("Count_")]
    X = df_bait[count_cols].to_numpy()
    return X, df_bait


def _summarize_histories_for_csv(histories):
    """
    Summarize EM histories for writing to a CSV file.

    Each row corresponds to one EM iteration and includes:
    log likelihood, mean lambda1, mean lambda2,
    mixture proportions pi1 and pi2,
    and mean responsibilities gamma1 and gamma2.
    """
    n_iter = len(histories["loglik_history"])
    rows = []
    for i in range(n_iter):
        lambda1_mean = histories["lambda1_history"][i].mean()
        lambda2_mean = histories["lambda2_history"][i].mean()
        pi_vec = histories["pi_history"][i]
        gamma_mean = histories["gamma_history"][i].mean(axis=0)
        row = {
            "iteration": i,
            "loglik": histories["loglik_history"][i],
            "lambda1_mean": lambda1_mean,
            "lambda2_mean": lambda2_mean,
            "pi1": pi_vec[0],
            "pi2": pi_vec[1],
            "gamma1_mean": gamma_mean[0],
            "gamma2_mean": gamma_mean[1]
        }
        rows.append(row)
    return pd.DataFrame(rows)


# %% Pipeline

def run_classical_pipeline(
    input_data,
    bait_names,
    metadata,
    max_iter=100,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=None,
    verbose=False,
    make_plots=True,
    plot_dir=None,
    mode="ide"
):
    """
    Run classical SAINT for multiple baits.

    lambda1 is the background Poisson rate for each prey (protein).
    lambda2 is the signal Poisson rate for each prey (protein).
    pi is the vector of mixture proportions for the two components.
    gamma represents the posterior probability that each prey belongs
    to each component.

    In mode, "ide", the function returns a dictionary keyed by bait name.
    Each entry contains:
        output_df: final lambda1, lambda2, pi, and gamma values
        merged_df: input data merged with output_df
        raw_outputs: histories of lambda, pi, gamma, and log likelihood
        metadata: bait and biological bait names

    In mode, "cli", the same structure is returned, but plots are not shown.
    Plots are saved to plot_dir if provided. For each bait, a bait-specific
    output .csv and histories .csv are written.

    If the input is in wide format (Protein by bait columns), the function
    performs implicit conversion to long format using extract_bait_matrix.
    """

    # Implicit conversion from wide format to long format if needed
    if "Bait" not in input_data.columns and "Protein" in input_data.columns:
        input_data = extract_bait_matrix(input_data)

    results = {}

    for bait in bait_names:

        biological_bait = metadata["bait_to_biological_name"][bait]

        df_bait = input_data[input_data["Bait"] == bait].copy()
        X, df_bait = reshape_counts_df(
            df_bait,
            bait,
            biological_bait
        )

        lambda1_init = np.maximum(X.mean(axis=1) * 0.5, 1e-3)
        lambda2_init = np.maximum(X.mean(axis=1) * 2.0, 1e-3)
        pi_init = np.array([0.8, 0.2])
        alpha = np.array([2.0, 2.0])

        hyperparams = {
            "lambda1_init": lambda1_init,
            "lambda2_init": lambda2_init,
            "pi_init": pi_init,
            "alpha": alpha
        }

        results_em = run_em(
            X,
            hyperparams,
            biological_bait,
            max_iter=max_iter,
            tol_loglik=tol_loglik,
            tol_params=tol_params,
            seed=seed,
            verbose=verbose
        )

        histories = {
            "loglik_history": results_em["loglik_history"],
            "lambda1_history": results_em["lambda1_history"],
            "lambda2_history": results_em["lambda2_history"],
            "pi_history": results_em["pi_history"],
            "gamma_history": results_em["gamma_history"]
        }

        if make_plots:
            figs = make_classical_plots(histories, bait)
            if mode == "ide":
                for fig in figs.values():
                    fig.show()
            if mode == "cli" and plot_dir is not None:
                os.makedirs(plot_dir, exist_ok=True)
                for name, fig in figs.items():
                    fig_path = os.path.join(
                        plot_dir,
                        f"{bait}_{name}.png"
                    )
                    fig.savefig(fig_path)
                plt.close("all")

        final_lambda1 = histories["lambda1_history"][-1]
        final_lambda2 = histories["lambda2_history"][-1]
        final_pi = histories["pi_history"][-1]
        final_gamma = histories["gamma_history"][-1]

        output_df = pd.DataFrame({
            "Protein": df_bait["Protein"].values,
            "lambda1": final_lambda1,
            "lambda2": final_lambda2,
            "pi1": final_pi[0],
            "pi2": final_pi[1],
            "gamma1": final_gamma[:, 0],
            "gamma2": final_gamma[:, 1]
        })

        merged_df = pd.merge(
            df_bait,
            output_df,
            on="Protein",
            how="left"
        )

        if mode == "cli":
            if plot_dir is not None:
                os.makedirs(plot_dir, exist_ok=True)
                out_path = os.path.join(plot_dir, f"{bait}_output.csv")
                output_df.to_csv(out_path, index=False)
                hist_df = _summarize_histories_for_csv(histories)
                hist_path = os.path.join(plot_dir, f"{bait}_histories.csv")
                hist_df.to_csv(hist_path, index=False)

        results[bait] = {
            "output_df": output_df,
            "merged_df": merged_df,
            "raw_outputs": histories,
            "metadata": {
                "bait": bait,
                "biological_bait": biological_bait
            }
        }

    return results









