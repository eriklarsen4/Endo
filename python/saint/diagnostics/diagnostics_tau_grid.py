# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 15:30:50 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% Summaries
def summarize_tau_grid(tau_grid_results):
    """
    Summarize the results of the tau grid search for a single bait.
    Returns a DataFrame with one row per tau value.
    """
    rows = []

    for tau, res in tau_grid_results.items():
        lambda1 = res["lambda1"]
        lambda2 = res["lambda2"]
        lambda3 = res["lambda3"]

        pi = res["pi"]              # shape (3,)
        gamma = res["gamma"]        # shape (n, 3)

        rows.append({
            "tau": tau,
            "loglik": res["loglik_history"][-1],
            "lambda1_mean": float(np.mean(lambda1)),
            "lambda2_mean": float(np.mean(lambda2)),
            "lambda3_mean": float(np.mean(lambda3)),
            "pi1": float(pi[0]),
            "pi2": float(pi[1]),
            "pi3": float(pi[2]),
            "mean_gamma3": float(np.mean(gamma[:, 2])),
        })

    return pd.DataFrame(rows).sort_values("tau")

# %% Plotting

def plot_tau_grid(df, best_tau, bait_name):
    """
    Create diagnostic plots for the tau grid search for a single bait.
    Returns a matplotlib Figure object.
    """

    sns.set_style("whitegrid")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Log-likelihood
    ax = axes[0, 0]
    ax.plot(df["tau"], df["loglik"], marker="o")
    ax.axvline(best_tau, color="red")
    ax.set_xlabel("Tau")
    ax.set_ylabel("Final log likelihood")
    ax.set_title(f"Tau vs Final Log-Likelihood — {bait_name}")

    # Lambda parameters (means)
    ax = axes[0, 1]
    ax.plot(df["tau"], df["lambda1_mean"], label="lambda1 mean")
    ax.plot(df["tau"], df["lambda2_mean"], label="lambda2 mean")
    ax.plot(df["tau"], df["lambda3_mean"], label="lambda3 mean")
    ax.axvline(best_tau, color="red")
    ax.set_xlabel("Tau")
    ax.set_ylabel("Lambda value")
    ax.set_title(f"Tau vs Lambda Parameters — {bait_name}")
    ax.legend()

    # Mixture weights
    ax = axes[1, 0]
    ax.plot(df["tau"], df["pi1"], label="pi1")
    ax.plot(df["tau"], df["pi2"], label="pi2")
    ax.plot(df["tau"], df["pi3"], label="pi3")
    ax.axvline(best_tau, color="red")
    ax.set_xlabel("Tau")
    ax.set_ylabel("Mixture weight")
    ax.set_title(f"Tau vs Mixture Weights — {bait_name}")
    ax.legend()

    # Mean gamma3 (signal occupancy)
    ax = axes[1, 1]
    ax.plot(df["tau"], df["mean_gamma3"], marker="o")
    ax.axvline(best_tau, color="red")
    ax.set_xlabel("Tau")
    ax.set_ylabel("Mean gamma3")
    ax.set_title(f"Tau vs Mean Gamma3 (Signal Occupancy) — {bait_name}")

    plt.tight_layout()
    return fig


# %% High-level entry point

def diagnostics_tau_grid(tau_info, bait_name):
    """
    Produce tau grid diagnostics for a single bait.
    Returns a dictionary containing:
        - summary_df: DataFrame of tau grid results
        - figure: matplotlib Figure object
        - best_tau: the selected tau value
    """

    tau_grid_results = tau_info["tau_grid_results"]
    best_tau = tau_info["best_tau"]

    df = summarize_tau_grid(tau_grid_results)
    fig = plot_tau_grid(df, best_tau, bait_name)

    return {
        "bait": bait_name,
        "best_tau": best_tau,
        "summary_df": df,
        "figure": fig
    }