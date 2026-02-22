# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 18:52:18 2026

@author: Erik
"""

# %% Imports

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# %% Classical diagnostics

def make_classical_plots(results_em, bait_name):
    """
    Create diagnostic plots for the classical SAINT EM results.
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
    plt.title(f"Classical SAINT log likelihood trajectory for {bait_name}")
    figs["loglik"] = fig_loglik

    # %% Lambda trajectories
    fig_lambda = plt.figure()
    lambda1_hist = np.array(results_em["lambda1_history"])
    lambda2_hist = np.array(results_em["lambda2_history"])

    plt.plot(lambda1_hist.mean(axis=1), label="lambda1 mean")
    plt.plot(lambda2_hist.mean(axis=1), label="lambda2 mean")
    plt.xlabel("Iteration")
    plt.ylabel("Mean lambda value")
    plt.title(f"Classical SAINT lambda trajectories for {bait_name}")
    plt.legend()
    figs["lambda"] = fig_lambda

    # %% Pi trajectory
    fig_pi = plt.figure()
    pi_hist = np.array(results_em["pi_history"])
    plt.plot(pi_hist[:, 0], label="pi1")
    plt.plot(pi_hist[:, 1], label="pi2")
    plt.xlabel("Iteration")
    plt.ylabel("Mixture weight")
    plt.title(f"Classical SAINT pi trajectories for {bait_name}")
    plt.legend()
    figs["pi"] = fig_pi

    # %% Gamma distributions
    fig_gamma = plt.figure()
    gamma_final = results_em["gamma"]
    plt.hist(gamma_final[:, 0], bins=30, alpha=0.5, label="gamma1")
    plt.hist(gamma_final[:, 1], bins=30, alpha=0.5, label="gamma2")
    plt.xlabel("Gamma value")
    plt.ylabel("Count")
    plt.title(f"Classical SAINT gamma distributions for {bait_name}")
    plt.legend()
    figs["gamma"] = fig_gamma

    return figs

