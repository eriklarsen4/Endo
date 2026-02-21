# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 13:17:16 2026

@author: Erik
"""

# %% Imports

import numpy as np
import matplotlib.pyplot as plt

# %% Parameter trajectory plots

def plot_lambda_trajectories(histories):
    """
    Plot lambda trajectories for each component.
    
    Parameters
    histories is a dictionary returned by the EM loop
    
    Returns
    A figure with lambda trajectories
    """

    lambda1_hist = histories["lambda1"]
    lambda2_hist = histories["lambda2"]
    lambda3_hist = histories["lambda3"]

    fig, ax = plt.subplots()
    ax.plot([np.mean(v) for v in lambda1_hist], label="lambda one")
    ax.plot([np.mean(v) for v in lambda2_hist], label="lambda two")
    ax.plot([np.mean(v) for v in lambda3_hist], label="lambda three")
    ax.set_xlabel("iteration")
    ax.set_ylabel("mean lambda")
    ax.legend()
    return fig


def plot_pi_trajectory(histories):
    """
    Plot pi trajectory.
    
    Parameters
    histories is a dictionary returned by the EM loop
    
    Returns
    A figure with pi trajectories
    """

    pi_hist = histories["pi"]

    fig, ax = plt.subplots()
    ax.plot([v[0] for v in pi_hist], label="pi one")
    ax.plot([v[1] for v in pi_hist], label="pi two")
    ax.plot([v[2] for v in pi_hist], label="pi three")
    ax.set_xlabel("iteration")
    ax.set_ylabel("pi value")
    ax.legend()
    return fig


def plot_tau_trajectory(histories):
    """
    Plot tau trajectory.
    
    Parameters
    histories is a dictionary returned by the EM loop
    
    Returns
    A figure with tau trajectory
    """

    tau_hist = histories["tau"]

    fig, ax = plt.subplots()
    ax.plot(tau_hist, label="tau")
    ax.set_xlabel("iteration")
    ax.set_ylabel("tau value")
    ax.legend()
    return fig


def build_parameter_trajectories(histories):
    """
    Build all parameter trajectory figures.
    
    Parameters
    histories is a dictionary returned by the EM loop
    
    Returns
    A dictionary of figures
    """

    figs = {}
    figs["lambda"] = plot_lambda_trajectories(histories)
    figs["pi"] = plot_pi_trajectory(histories)
    figs["tau"] = plot_tau_trajectory(histories)

    return figs