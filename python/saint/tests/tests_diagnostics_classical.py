# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:59:59 2026

@author: Erik
"""

# %% Imports

import numpy as np
from matplotlib.figure import Figure

from saint.diagnostics.diagnostics_classical import make_classical_plots


# %% Fixtures

def make_raw_outputs():
    """
    Create a minimal raw output dictionary for testing classical diagnostics.
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
        "column_labels": ["C1"],
        "mapping": {"Protein": [0, 1, 2]}
    }


# %% Tests

def test_classical_diagnostics_return_figures():
    """
    Test that classical diagnostics return a dictionary of figures.
    """

    raw = make_raw_outputs()
    figs = make_classical_plots(raw, output_dir=None)

    assert isinstance(figs, dict)
    assert len(figs) > 0

    for f in figs.values():
        assert isinstance(f, Figure)


def test_classical_diagnostics_do_not_modify_raw_outputs():
    """
    Test that classical diagnostics do not modify the raw output dictionary.
    """

    raw = make_raw_outputs()
    before = raw["lambda1"].copy()

    _ = make_classical_plots(raw, output_dir=None)

    after = raw["lambda1"]
    assert np.allclose(before, after)