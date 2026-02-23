# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:59 2026

@author: Erik
"""

# %% Imports
import numpy as np
from saint.diagnostics.diagnostics_classical import make_classical_plots
import matplotlib.figure as mpl_fig


# %% Tests

def test_make_classical_plots_returns_figures():
    # Minimal synthetic EM results
    results_em = {
        "loglik_history": [0.0, 1.0, 2.0],
        "lambda1_history": [np.array([1.0, 2.0]), np.array([1.5, 2.5])],
        "lambda2_history": [np.array([3.0, 4.0]), np.array([3.5, 4.5])],
        "pi_history": [np.array([0.7, 0.3]), np.array([0.6, 0.4])],
        "gamma_history": [np.array([[0.8, 0.2], [0.3, 0.7]])],
        "gamma": np.array([[0.8, 0.2], [0.3, 0.7]]),
    }

    figs = make_classical_plots(results_em, bait_name="B1")

    # Expected figure keys
    expected_keys = {"loglik", "lambda", "pi", "gamma"}
    assert expected_keys.issubset(figs.keys())

    # All returned objects must be matplotlib Figure instances
    for fig in figs.values():
        assert isinstance(fig, mpl_fig.Figure)