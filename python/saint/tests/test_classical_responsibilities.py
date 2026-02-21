# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:28:38 2026

@author: Erik
"""

# %% Imports

import numpy as np
from saint.model.classical_responsibilities import compute_responsibilities


# %% Tests

def test_gamma_shape_two_components():
    X = np.array([[5, 3], [0, 1]])
    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    pi = np.array([0.7, 0.3])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    assert gamma.shape == (2, 2)


def test_gamma_rows_sum_to_one():
    X = np.array([[4, 4], [2, 2], [0, 1]])
    lambda1 = np.array([1.0, 1.0, 1.0])
    lambda2 = np.array([3.0, 3.0, 3.0])
    pi = np.array([0.5, 0.5])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    row_sums = gamma.sum(axis=1)
    assert np.allclose(row_sums, 1.0)


def test_gamma_symmetry_when_rates_equal():
    X = np.array([[5, 5]])
    lambda1 = np.array([2.0])
    lambda2 = np.array([2.0])
    pi = np.array([0.5, 0.5])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    assert np.allclose(gamma[0, 0], gamma[0, 1])


def test_gamma_extreme_rate_difference():
    X = np.array([[10, 10]])
    lambda1 = np.array([1.0])
    lambda2 = np.array([20.0])
    pi = np.array([0.5, 0.5])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    assert gamma[0, 1] > gamma[0, 0]


def test_no_nans_or_infs():
    X = np.array([[3, 1], [0, 0]])
    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    pi = np.array([0.6, 0.4])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    assert np.isfinite(gamma).all()