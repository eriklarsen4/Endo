# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 15:32:08 2026

@author: Erik
"""

# %% Imports

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde


# %% Gamma3 density across baits

def plot_gamma3_density_across_baits(results_by_bait, positive_control_bait=None):
    """
    Plot the density of mean gamma3 values across all baits.
    Optionally highlight a positive control bait.
    
    Parameters
    
    results_by_bait : dict
        Dictionary mapping bait_name -> EM result dict.
        Each EM result dict must contain "gamma" (n × 3 array).
    positive_control_bait : str or None
        Name of a bait to highlight in the density plot.
    
    Returns
    
    fig : matplotlib Figure
        The density plot figure of the average posterior probability of the signal component (gamma_3).
    """

    sns.set_style("whitegrid")

    # Extract mean gamma3 for each bait
    bait_names = []
    mean_gamma3_values = []

    for bait, res in results_by_bait.items():
        gamma = res["gamma"]  # n × 3
        mean_gamma3 = gamma[:, 2].mean()  # gamma3
        bait_names.append(bait)
        mean_gamma3_values.append(mean_gamma3)

    mean_gamma3_values = np.array(mean_gamma3_values)

    # KDE density estimate
    kde = gaussian_kde(mean_gamma3_values)
    xs = np.linspace(0, 1, 400)
    ys = kde(xs)

    # Plot
    fig = plt.figure(figsize=(8, 5))
    plt.plot(xs, ys, label="Density of mean gamma3 across baits")
    plt.xlabel("Mean gamma3")
    plt.ylabel("Density")
    plt.title("Distribution of Mean Gamma3 Across Baits")

    # Highlight positive control bait if provided
    if positive_control_bait is not None:
        pc_value = results_by_bait[positive_control_bait]["gamma"][:, 2].mean()
        plt.axvline(pc_value, color="red", linewidth=2,
                    label=f"{positive_control_bait} (positive control)")

    plt.legend()
    return fig


# %% Summary table across baits

def summarize_gamma3_across_baits(results_by_bait):
    """
    Return a simple summary table of mean gamma3 values across baits.
    
    Parameters
    
    results_by_bait : dict
        Dictionary mapping bait_name -> EM result dict.
    
    Returns
    
    summary : dict
        Dictionary mapping bait_name -> mean_gamma3.
    """

    summary = {}

    for bait, res in results_by_bait.items():
        gamma = res["gamma"]
        summary[bait] = gamma[:, 2].mean()  # mean gamma3

    return summary


# %% High-level entry point

def diagnostics_pipeline(results_by_bait, positive_control_bait=None):
    """
    Produce pipeline-level diagnostics across all baits.
    Includes:
        - density plot of mean gamma3 across baits
        - summary table of mean gamma3 values
    
    Parameters
    
    results_by_bait : dict
        Dictionary mapping bait_name -> EM result dict.
    positive_control_bait : str or None
        Name of a bait to highlight in the density plot.
    
    Returns
    
    dict containing:
        - gamma3_density_figure
        - gamma3_summary
    """

    fig_density = plot_gamma3_density_across_baits(
        results_by_bait,
        positive_control_bait=positive_control_bait
    )

    summary = summarize_gamma3_across_baits(results_by_bait)

    return {
        "gamma3_density_figure": fig_density,
        "gamma3_summary": summary
    }