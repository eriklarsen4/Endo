# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 14:51:03 2026

@author: Erik
"""

# %% Imports

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import ScalarFormatter


# %% Hierarchical Plots
def make_hierarchical_plots(results_em, bait_name):
    """
    Faceted hierarchical diagnostics: trajectories in a 2×2 grid.
    Returns a dict of figures.

    This function expects the EM result dictionary produced by
    run_em_hierarchical (or extracted from the tau-grid / pipeline),
    with the following keys:

        loglik_history : list of float
        lambda1_history, lambda2_history, lambda3_history : list of arrays
        pi_history : list of length-3 arrays
        alpha_history : list of length-3 arrays
        a_history, b_history : list of length-3 arrays
        gamma : final n × 3 responsibility matrix
    """
    
    sns.set_style("whitegrid")
    figs = {}

    # Extract histories
    loglik = np.array(results_em["loglik_history"])
    lambda1 = np.array(results_em["lambda1_history"]).mean(axis=1)
    lambda2 = np.array(results_em["lambda2_history"]).mean(axis=1)
    lambda3 = np.array(results_em["lambda3_history"]).mean(axis=1)
    pi_hist = np.array(results_em["pi_history"])
    alpha_hist = np.array(results_em["alpha_history"])
    a_hist = np.array(results_em["a_history"])
    b_hist = np.array(results_em["b_history"])

    # --- Faceted trajectories ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Log-likelihood
    ax = axes[0, 0]
    ax.plot(loglik, marker="o")
    ax.set_title(f"Log-likelihood — {bait_name}")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Log-likelihood")
    ax.ticklabel_format(style="plain", axis="both")

    # Lambda means
    ax = axes[0, 1]
    ax.plot(lambda1, label="λ1 mean")
    ax.plot(lambda2, label="λ2 mean")
    ax.plot(lambda3, label="λ3 mean")
    ax.set_title(f"Lambda trajectories — {bait_name}")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Mean λ")
    ax.ticklabel_format(style="plain", axis="both")
    ax.legend()

    # Mixture weights
    ax = axes[1, 0]
    ax.plot(pi_hist[:, 0], label="π1")
    ax.plot(pi_hist[:, 1], label="π2")
    ax.plot(pi_hist[:, 2], label="π3")
    ax.set_title(f"Mixture weights — {bait_name}")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("π")
    ax.ticklabel_format(style="plain", axis="both")
    ax.legend()

    # Prior parameters
    ax = axes[1, 1]
    ax.plot(alpha_hist[:, 0], label="α1")
    ax.plot(alpha_hist[:, 1], label="α2")
    ax.plot(alpha_hist[:, 2], label="α3")
    ax.plot(a_hist[:, 0], label="a1")
    ax.plot(a_hist[:, 1], label="a2")
    ax.plot(a_hist[:, 2], label="a3")
    ax.plot(b_hist[:, 0], label="b1")
    ax.plot(b_hist[:, 1], label="b2")
    ax.plot(b_hist[:, 2], label="b3")
    ax.set_title(f"Prior parameters — {bait_name}")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Value")
    ax.ticklabel_format(style="plain", axis="both")
    ax.legend()

    plt.tight_layout()
    figs["trajectories"] = fig
    
    # --- Gamma distributions ---
    fig_gamma = plt.figure(figsize=(6, 4))
    gamma_final = results_em["gamma"]
    plt.hist(gamma_final[:, 0], bins=30, alpha=0.5, label="γ1")
    plt.hist(gamma_final[:, 1], bins=30, alpha=0.5, label="γ2")
    plt.hist(gamma_final[:, 2], bins=30, alpha=0.5, label="γ3")
    plt.xlabel("Gamma value")
    plt.ylabel("Count")
    plt.title(f"Gamma distributions — {bait_name}")
    plt.legend()
    figs["gamma"] = fig_gamma

    return figs

# %% KDE of average gamma3 posterior probabilities

def plot_gamma3_density(results_df, protein_of_interest="TMEM184B", ax=None):
    """
    KDE of average gamma3 across TMEMV5 and TMEMmyc,
    with a vertical red line at the specified protein's average gamma3.

    Parameters
    
    results_df : DataFrame
        Output results_df from the hierarchical pipeline, with columns:
        Protein, bait, lambda1, lambda2, lambda3, tau,
        pi1, pi2, pi3, gamma1, gamma2, gamma3.

    protein_of_interest : str
        Protein identifier whose average gamma3 (across TMEMV5 and TMEMmyc)
        will be highlighted with a vertical line, if present.

    ax : matplotlib.axes.Axes, optional
        Existing axis to draw on. If None, a new figure and axis are created.

    Returns
    
    ax : matplotlib.axes.Axes
        Axis containing the KDE plot.
    """

    df = results_df  # pipeline passes results_df directly

    # Extract gamma3 for each treatment bait
    gamma3_v5 = df[df["bait"] == "TMEMV5"][["Protein", "gamma3"]].rename(
        columns={"gamma3": "gamma3_v5"}
    )
    gamma3_myc = df[df["bait"] == "TMEMmyc"][["Protein", "gamma3"]].rename(
        columns={"gamma3": "gamma3_myc"}
    )

    # Merge and compute average gamma3
    merged = gamma3_v5.merge(gamma3_myc, on="Protein", how="inner")
    merged["gamma3_avg"] = (merged["gamma3_v5"] + merged["gamma3_myc"]) / 2.0

    # Extract the value for the protein of interest
    if protein_of_interest in merged["Protein"].values:
        marker_value = float(
            merged.loc[merged["Protein"] == protein_of_interest, "gamma3_avg"]
        )
    else:
        marker_value = None  # no vertical line if protein not found

    # Create axis if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    # KDE plot
    vals = merged["gamma3_avg"].round(6)
    sns.kdeplot(vals, fill=True, color="steelblue", ax=ax)

    if marker_value is not None:
        marker_value = float(round(marker_value, 6))
        ax.axvline(marker_value, color="red", linewidth=2)
        
    fmt = ScalarFormatter(useOffset=False)
    fmt.set_powerlimits((-3, 4))
    fmt.set_useMathText(False)
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)


    ax.set_xlabel("Average gamma3 value")
    ax.set_ylabel("Density")
    ax.set_title("Combined gamma3 density for TMEMV5 and TMEMmyc")

    # Disable scientific notation
    ax.ticklabel_format(style="plain", axis="both")
    ax.get_xaxis().get_offset_text().set_visible(False)
    ax.get_yaxis().get_offset_text().set_visible(False)

    return ax


# %% Pipeline-level KDE of posterior probabilities of gamma3
def plot_multi_bait_summary(results_df, bait_list):
    """
    Plot mean gamma3 per bait as a simple multi-bait summary.

    Parameters
    
    results_df : DataFrame
        Output results_df from the hierarchical pipeline.

    bait_list : list of str
        Ordered list of bait names to include and display on the x-axis.
    """
    mean_gamma3 = results_df.groupby("bait")["gamma3"].mean().reindex(bait_list)

    plt.figure(figsize=(8, 4))
    plt.plot(bait_list, mean_gamma3, marker="o")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean gamma3")
    plt.title("Mean gamma3 per bait")
    plt.grid(True, alpha=0.3)