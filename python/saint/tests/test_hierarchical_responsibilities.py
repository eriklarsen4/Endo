# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:29:14 2026

@author: Erik
"""

# %% Imports

import numpy as np
from saint.model.hierarchical_responsibilities import compute_responsibilities


# %% Tests

def test_gamma_shape_three_components():
    X = np.array([[5, 3], [1, 0]])
    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    lambda3 = np.array([4.0, 4.0])
    pi = np.array([0.5, 0.3, 0.2])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    assert gamma.shape == (2, 3)


def test_gamma_rows_sum_to_one():
    X = np.array([[4, 4], [2, 2], [0, 1]])
    lambda1 = np.array([1.0, 1.0, 1.0])
    lambda2 = np.array([3.0, 3.0, 3.0])
    lambda3 = np.array([6.0, 6.0, 6.0])
    pi = np.array([0.4, 0.4, 0.2])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    row_sums = gamma.sum(axis=1)
    assert np.allclose(row_sums, 1.0)


def test_gamma_symmetry_when_all_rates_equal():
    X = np.array([[5, 5]])
    lambda1 = np.array([2.0])
    lambda2 = np.array([2.0])
    lambda3 = np.array([2.0])
    pi = np.array([1/3, 1/3, 1/3])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    assert np.allclose(gamma[0], np.array([1/3, 1/3, 1/3]))


def test_gamma_extreme_rate_difference():
    X = np.array([[10, 10]])
    lambda1 = np.array([1.0])
    lambda2 = np.array([5.0])
    lambda3 = np.array([20.0])
    pi = np.array([0.3, 0.3, 0.4])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    assert gamma[0, 2] > gamma[0, 1] > gamma[0, 0]


def test_no_nans_or_infs():
    X = np.array([[3, 1], [0, 0]])
    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    lambda3 = np.array([3.0, 3.0])
    pi = np.array([0.5, 0.3, 0.2])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    assert np.isfinite(gamma).all()