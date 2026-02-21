# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 20:09:11 2026

@author: Erik
"""

# %% Imports

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# %% Log likelihood plot

def _plot_loglik(histories, output_dir):
    loglik = histories["loglik"]

    fig = plt.figure(figsize=(6, 4))
    sns.lineplot(x=np.arange(len(loglik)), y=loglik)
    plt.xlabel("iteration")
    plt.ylabel("log likelihood")
    plt.title("hierarchical SAINT log likelihood")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "loglik.png"))
    plt.close()
    return fig

# %% Lambda mean trajectories

def _plot_lambda_means(histories, output_dir):
    lambda1_hist = np.array(histories["lambda1"])
    lambda2_hist = np.array(histories["lambda2"])
    lambda3_hist = np.array(histories["lambda3"])

    fig = plt.figure(figsize=(6, 4))
    sns.lineplot(x=np.arange(lambda1_hist.shape[0]), y=lambda1_hist.mean(axis=1), label="lambda1 mean")
    sns.lineplot(x=np.arange(lambda2_hist.shape[0]), y=lambda2_hist.mean(axis=1), label="lambda2 mean")
    sns.lineplot(x=np.arange(lambda3_hist.shape[0]), y=lambda3_hist.mean(axis=1), label="lambda3 mean")
    plt.xlabel("iteration")
    plt.ylabel("mean lambda")
    plt.title("lambda trajectories")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "lambda_means.png"))
    plt.close()
    return fig

# %% Tau trajectory

def _plot_tau(histories, output_dir):
    tau_hist = histories["tau"]

    fig = plt.figure(figsize=(6, 4))
    sns.lineplot(x=np.arange(len(tau_hist)), y=tau_hist)
    plt.xlabel("iteration")
    plt.ylabel("tau")
    plt.title("tau trajectory")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "tau.png"))
    plt.close()
    return fig

# %% Gamma mean trajectories

def _plot_gamma_means(histories, output_dir):
    gamma_hist = np.array(histories["gamma"])
    gamma_means = gamma_hist.mean(axis=1)

    fig = plt.figure(figsize=(6, 4))
    sns.lineplot(x=np.arange(gamma_means.shape[0]), y=gamma_means[:, 0], label="gamma1 mean")
    sns.lineplot(x=np.arange(gamma_means.shape[0]), y=gamma_means[:, 1], label="gamma2 mean")
    sns.lineplot(x=np.arange(gamma_means.shape[0]), y=gamma_means[:, 2], label="gamma3 mean")
    plt.xlabel("iteration")
    plt.ylabel("mean gamma")
    plt.title("gamma trajectories")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "gamma_means.png"))
    plt.close()
    return fig

# %% Pi trajectories

def _plot_pi(histories, output_dir):
    pi_hist = np.array(histories["pi"])

    fig = plt.figure(figsize=(6, 4))
    sns.lineplot(x=np.arange(pi_hist.shape[0]), y=pi_hist[:, 0], label="pi1")
    sns.lineplot(x=np.arange(pi_hist.shape[0]), y=pi_hist[:, 1], label="pi2")
    sns.lineplot(x=np.arange(pi_hist.shape[0]), y=pi_hist[:, 2], label="pi3")
    plt.xlabel("iteration")
    plt.ylabel("pi value")
    plt.title("pi trajectories")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "pi_trajectories.png"))
    plt.close()
    return fig

# %% Gamma distributions

def _plot_gamma_distributions(histories, output_dir):
    gamma_hist = histories["gamma"]
    final_gamma = gamma_hist[-1]

    fig = plt.figure(figsize=(6, 4))
    sns.histplot(final_gamma[:, 0], kde=True, label="gamma1", color="blue")
    sns.histplot(final_gamma[:, 1], kde=True, label="gamma2", color="orange")
    sns.histplot(final_gamma[:, 2], kde=True, label="gamma3", color="green")
    plt.xlabel("gamma value")
    plt.ylabel("count")
    plt.title("final gamma distributions")
    plt.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "gamma_distributions.png"))
    plt.close()
    return fig

# %% Tau distribution

def _plot_tau_distribution(histories, output_dir):
    tau_hist = histories["tau"]

    fig = plt.figure(figsize=(6, 4))
    sns.histplot(tau_hist, kde=True)
    plt.xlabel("tau")
    plt.ylabel("count")
    plt.title("tau distribution across iterations")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "tau_distribution.png"))
    plt.close()
    return fig

# %% Lambda distributions

def _plot_lambda_distributions(histories, output_dir):
    lambda1_hist = np.array(histories["lambda1"])
    lambda2_hist = np.array(histories["lambda2"])
    lambda3_hist = np.array(histories["lambda3"])

    final_lambda1 = lambda1_hist[-1]
    final_lambda2 = lambda2_hist[-1]
    final_lambda3 = lambda3_hist[-1]

    fig = plt.figure(figsize=(6, 4))
    sns.histplot(final_lambda1, kde=True, label="lambda1", color="blue")
    sns.histplot(final_lambda2, kde=True, label="lambda2", color="orange")
    sns.histplot(final_lambda3, kde=True, label="lambda3", color="green")
    plt.xlabel("lambda value")
    plt.ylabel("count")
    plt.title("final lambda distributions")
    plt.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "lambda_distributions.png"))
    plt.close()
    return fig

# %% Gamma per bait

def _plot_gamma_per_bait(raw_outputs, output_dir):
    histories = raw_outputs["histories"]
    mapping = raw_outputs["mapping"]
    gamma_final = histories["gamma"][-1]

    by_bait = mapping["by_bait"]
    column_labels = raw_outputs["column_labels"]

    fig = plt.figure(figsize=(7, 5))

    for bait, cols in by_bait.items():
        indices = [column_labels.index(c) for c in cols]
        gamma_subset = gamma_final[:, 0]
        sns.kdeplot(gamma_subset, label=bait)

    plt.xlabel("gamma1 value")
    plt.ylabel("density")
    plt.title("gamma1 distributions per bait")
    plt.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "gamma_per_bait.png"))
    plt.close()
    return fig

# %% Public entry point

def make_hierarchical_plots(raw_outputs, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    histories = raw_outputs["histories"]

    figs = {}
    figs["loglik"] = _plot_loglik(histories, output_dir)
    figs["lambda_means"] = _plot_lambda_means(histories, output_dir)
    figs["tau"] = _plot_tau(histories, output_dir)
    figs["gamma_means"] = _plot_gamma_means(histories, output_dir)
    figs["pi"] = _plot_pi(histories, output_dir)
    figs["gamma_distributions"] = _plot_gamma_distributions(histories, output_dir)
    figs["tau_distribution"] = _plot_tau_distribution(histories, output_dir)
    figs["lambda_distributions"] = _plot_lambda_distributions(histories, output_dir)
    figs["gamma_per_bait"] = _plot_gamma_per_bait(raw_outputs, output_dir)

    return figs