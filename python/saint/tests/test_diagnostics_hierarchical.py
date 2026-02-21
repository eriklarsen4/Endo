# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:19:02 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from saint.diagnostics.diagnostics_hierarchical import make_hierarchical_plots


# %% Fixtures

def make_raw_outputs():
    """
    Create a minimal raw output dictionary for testing diagnostics.
    """

    lambda1 = np.array([1.0, 2.0, 3.0])
    lambda2 = np.array([0.5, 0.5, 0.5])
    lambda3 = np.array([0.1, 0.1, 0.1])

    pi = np.array([0.33, 0.33, 0.34])

    gamma = np.array([
        [0.8, 0.1, 0.1],
        [0.2, 0.3, 0.5],
        [0.1, 0.1, 0.8]
    ])

    tau = np.array([1.0, 1.0, 1.0])

    histories = {
        "lambda1": [lambda1],
        "lambda2": [lambda2],
        "lambda3": [lambda3],
        "pi": [pi],
        "loglik": [-10.0, -5.0, -3.0]
    }

    return {
        "lambda1": lambda1,
        "lambda2": lambda2,
        "lambda3": lambda3,
        "pi": pi,
        "gamma": gamma,
        "tau": tau,
        "histories": histories,
        "column_labels": ["C1", "C2", "C3"],
        "mapping": {"Protein": [0, 1, 2]}
    }


# %% Tests

def test_make_plots_returns_figures():
    """
    Test that diagnostics return a dictionary of figures.
    """

    raw_outputs = make_raw_outputs()
    figs = make_hierarchical_plots(raw_outputs, output_dir=None)

    assert isinstance(figs, dict)
    assert len(figs) > 0

    for f in figs.values():
        assert isinstance(f, Figure)


def test_diagnostics_do_not_modify_raw_outputs():
    """
    Test that diagnostics do not change the raw output dictionary.
    """

    raw_outputs = make_raw_outputs()
    before = raw_outputs["lambda1"].copy()

    _ = make_hierarchical_plots(raw_outputs, output_dir=None)

    after = raw_outputs["lambda1"]
    assert np.allclose(before, after)