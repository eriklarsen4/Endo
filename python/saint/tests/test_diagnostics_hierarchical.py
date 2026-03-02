# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:19:02 2026

@author: Erik
"""

# %% Imports
import numpy as np
from saint.diagnostics.diagnostics_hierarchical import make_hierarchical_plots
import matplotlib.figure as mpl_fig


# %% Tests

def test_make_hierarchical_plots_returns_figures():
    # Minimal synthetic hierarchical EM results
    results_em = {
        "loglik_history": [0.0, 1.0, 2.0],

        "lambda1_history": [np.array([1.0, 2.0]), np.array([1.5, 2.5])],
        "lambda2_history": [np.array([3.0, 4.0]), np.array([3.5, 4.5])],
        "lambda3_history": [np.array([5.0, 6.0]), np.array([5.5, 6.5])],

        "pi_history": [
            np.array([0.6, 0.3, 0.1]),
            np.array([0.5, 0.3, 0.2]),
        ],

        "alpha_history": [
            np.array([2.0, 2.0, 2.0]),
            np.array([2.1, 2.0, 1.9]),
        ],
        "a_history": [
            np.array([1.0, 1.0, 1.0]),
            np.array([1.1, 1.0, 0.9]),
        ],
        "b_history": [
            np.array([0.5, 0.5, 0.5]),
            np.array([0.6, 0.5, 0.4]),
        ],

        "gamma_history": [
            np.array([[0.8, 0.1, 0.1], [0.2, 0.3, 0.5]])
        ],

        # Final gamma matrix for histogram
        "gamma": np.array([
            [0.8, 0.1, 0.1],
            [0.2, 0.3, 0.5],
        ]),
    }

    figs = make_hierarchical_plots(results_em, bait_name="B1")

    expected_keys = {
        "loglik",
        "lambda",
        "pi",
        "priors",
        "gamma",
        "gamma_heatmap",
    }
    assert expected_keys.issubset(figs.keys())

    # All returned objects must be matplotlib Figure instances
    for fig in figs.values():
        assert isinstance(fig, mpl_fig.Figure)