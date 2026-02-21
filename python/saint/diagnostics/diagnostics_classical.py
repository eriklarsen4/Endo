# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 21:45:43 2026

@author: Erik
"""

# %% Imports

import matplotlib.pyplot as plt
import seaborn as sns


# %% Plot helpers

def _plot_loglik(history, bait):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots()
    ax.plot(history, marker="o")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Log likelihood")
    ax.set_title(f"Log likelihood progression for bait {bait}")
    return fig


def _plot_lambda(history, label, bait):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots()
    for i, vec in enumerate(history):
        ax.plot(vec, label=f"iter {i}")
    ax.set_xlabel("Prey index")
    ax.set_ylabel(label)
    ax.set_title(f"{label} progression for bait {bait}")
    return fig


def _plot_pi(history, bait):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots()
    pi1 = [vec[0] for vec in history]
    pi2 = [vec[1] for vec in history]
    ax.plot(pi1, label="pi1", marker="o")
    ax.plot(pi2, label="pi2", marker="o")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("pi values")
    ax.set_title(f"pi progression for bait {bait}")
    ax.legend()
    return fig


def _plot_gamma(history, bait):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots()
    for i, mat in enumerate(history):
        ax.plot(mat.mean(axis=0), label=f"iter {i}")
    ax.set_xlabel("Component index")
    ax.set_ylabel("Mean gamma")
    ax.set_title(f"Gamma progression for bait {bait}")
    return fig


# %% Public entry point

def make_classical_plots(histories, bait):
    """
    Create diagnostic plots for classical SAINT.

    histories is a dictionary containing:
    loglik_history, lambda1_history, lambda2_history,
    pi_history, and gamma_history.

    Returns a dictionary of matplotlib Figure objects.
    """

    figs = {}

    figs["loglik"] = _plot_loglik(histories["loglik_history"], bait)
    figs["lambda1"] = _plot_lambda(histories["lambda1_history"], "lambda1", bait)
    figs["lambda2"] = _plot_lambda(histories["lambda2_history"], "lambda2", bait)
    figs["pi"] = _plot_pi(histories["pi_history"], bait)
    figs["gamma"] = _plot_gamma(histories["gamma_history"], bait)

    return figs


