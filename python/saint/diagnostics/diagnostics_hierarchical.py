# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 18:51:52 2026

@author: Erik
"""

# %% Imports

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# %% Hierarchical diagnostics

def make_hierarchical_plots(results_em, bait_name):
    """
    Create diagnostic plots for the hierarchical SAINT EM results.
    This function returns a dictionary of matplotlib Figure objects.
    Figures will automatically appear in IDEs such as Spyder or VS Code.
    """

    sns.set_style("whitegrid")

    figs = {}

    # %% Log likelihood trajectory
    fig_loglik = plt.figure()
    plt.plot(results_em["loglik_history"], marker="o")
    plt.xlabel("Iteration")
    plt.ylabel("Log likelihood")
    plt.title(f"Hierarchical SAINT log likelihood trajectory for {bait_name}")
    figs["loglik"] = fig_loglik

    # %% Lambda trajectories
    fig_lambda = plt.figure()
    lambda1_hist = np.array(results_em["lambda1_history"])
    lambda2_hist = np.array(results_em["lambda2_history"])
    lambda3_hist = np.array(results_em["lambda3_history"])

    plt.plot(lambda1_hist.mean(axis=1), label="lambda1 mean")
    plt.plot(lambda2_hist.mean(axis=1), label="lambda2 mean")
    plt.plot(lambda3_hist.mean(axis=1), label="lambda3 mean")
    plt.xlabel("Iteration")
    plt.ylabel("Mean lambda value")
    plt.title(f"Hierarchical SAINT lambda trajectories for {bait_name}")
    plt.legend()
    figs["lambda"] = fig_lambda

    # %% Pi trajectory
    fig_pi = plt.figure()
    pi_hist = np.array(results_em["pi_history"])
    plt.plot(pi_hist[:, 0], label="pi1")
    plt.plot(pi_hist[:, 1], label="pi2")
    plt.plot(pi_hist[:, 2], label="pi3")
    plt.xlabel("Iteration")
    plt.ylabel("Mixture weight")
    plt.title(f"Hierarchical SAINT pi trajectories for {bait_name}")
    plt.legend()
    figs["pi"] = fig_pi

    # %% Prior trajectories alpha, a, b
    fig_priors = plt.figure()
    alpha_hist = np.array(results_em["alpha_history"])
    a_hist = np.array(results_em["a_history"])
    b_hist = np.array(results_em["b_history"])

    plt.plot(alpha_hist[:, 0], label="alpha1")
    plt.plot(alpha_hist[:, 1], label="alpha2")
    plt.plot(alpha_hist[:, 2], label="alpha3")
    plt.plot(a_hist[:, 0], label="a1")
    plt.plot(a_hist[:, 1], label="a2")
    plt.plot(a_hist[:, 2], label="a3")
    plt.plot(b_hist[:, 0], label="b1")
    plt.plot(b_hist[:, 1], label="b2")
    plt.plot(b_hist[:, 2], label="b3")

    plt.xlabel("Iteration")
    plt.ylabel("Prior parameter value")
    plt.title(f"Hierarchical SAINT prior trajectories for {bait_name}")
    plt.legend()
    figs["priors"] = fig_priors

    # %% Gamma distributions
    fig_gamma = plt.figure()
    gamma_final = results_em["gamma"]
    plt.hist(gamma_final[:, 0], bins=30, alpha=0.5, label="gamma1")
    plt.hist(gamma_final[:, 1], bins=30, alpha=0.5, label="gamma2")
    plt.hist(gamma_final[:, 2], bins=30, alpha=0.5, label="gamma3")
    plt.xlabel("Gamma value")
    plt.ylabel("Count")
    plt.title(f"Hierarchical SAINT gamma distributions for {bait_name}")
    plt.legend()
    figs["gamma"] = fig_gamma

    return figs

# %% Distribution of posterior probabilities across baits
def plot_gamma3_density(results_df):
    """
    Plot the distribution of gamma3 values across all baits.
    """
    plt.figure()
    plt.hist(results_df["gamma3"], bins=50, alpha=0.7, color="steelblue")
    plt.xlabel("gamma3")
    plt.ylabel("Count")
    plt.title("Distribution of gamma3 across all baits")
    plt.grid(True, alpha=0.3)

# %% Pipeline-level kde of posterior probabililties of gamma3's
def plot_multi_bait_summary(results_df, bait_list):
    """
    Plot mean gamma3 per bait as a simple multi-bait summary.
    """
    mean_gamma3 = results_df.groupby("bait")["gamma3"].mean().reindex(bait_list)

    plt.figure(figsize=(8, 4))
    plt.plot(bait_list, mean_gamma3, marker="o")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Mean gamma3")
    plt.title("Mean gamma3 per bait")
    plt.grid(True, alpha=0.3)

