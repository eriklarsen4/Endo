# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:20:41 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pytest

from saint.validation.validate_hyperparams import validate_hyperparams


# %% Tests

def test_validate_hyperparams_accepts_valid_input():
    """
    Test that valid hyperparameters pass validation.
    """

    lambda1 = np.array([1.0, 2.0, 3.0])
    lambda2 = np.array([0.5, 0.5, 0.5])
    lambda3 = np.array([0.1, 0.1, 0.1])

    pi = np.array([0.33, 0.33, 0.34])

    validate_hyperparams(
        max_iter=10,
        tol=1e-4,
        lambda_vecs=[lambda1, lambda2, lambda3],
        pi=pi
    )


def test_validate_hyperparams_rejects_invalid_pi():
    """
    Test that invalid pi vectors raise an error.
    """

    lambda1 = np.array([1.0, 2.0, 3.0])
    lambda2 = np.array([0.5, 0.5, 0.5])
    lambda3 = np.array([0.1, 0.1, 0.1])

    pi = np.array([0.5, 0.5, 0.5])

    with pytest.raises(Exception):
        validate_hyperparams(
            max_iter=10,
            tol=1e-4,
            lambda_vecs=[lambda1, lambda2, lambda3],
            pi=pi
        )


def test_validate_hyperparams_rejects_mismatched_lambda_shapes():
    """
    Test that mismatched lambda vector shapes raise an error.
    """

    lambda1 = np.array([1.0, 2.0])
    lambda2 = np.array([0.5, 0.5, 0.5])
    lambda3 = np.array([0.1, 0.1, 0.1])

    pi = np.array([0.33, 0.33, 0.34])

    with pytest.raises(Exception):
        validate_hyperparams(
            max_iter=10,
            tol=1e-4,
            lambda_vecs=[lambda1, lambda2, lambda3],
            pi=pi
        )