# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 13:55:20 2026

@author: Erik
"""

# %% Imports
import numpy as np
from scipy.special import gammaln

from saint.model.classical_likelihood import compute_loglik


# %% Tests

def test_compute_classical_loglik_simple_case():
    """
    Minimal test verifying that compute_loglik matches a manually
    computed mixture-of-Poissons likelihood on a tiny dataset.
    """

    # X: 2 preys × 2 replicates
    X = np.array([
        [1, 2],
        [0, 1],
    ], dtype=float)

    # Classical SAINT uses only two components
    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])

    # Mixture weights
    pi = np.array([0.6, 0.4])

    # ----- Manual expected log-likelihood -----
    expected = 0.0
    for p in range(X.shape[0]):
        Xp = X[p]

        ll1 = (Xp * np.log(lambda1[p]) - lambda1[p] - gammaln(Xp + 1)).sum()
        ll2 = (Xp * np.log(lambda2[p]) - lambda2[p] - gammaln(Xp + 1)).sum()

        joint = pi[0] * np.exp(ll1) + pi[1] * np.exp(ll2)
        expected += np.log(joint)

    # ----- Actual -----
    loglik = compute_loglik(
        X,
        lambda1,
        lambda2,
        pi,
    )

    assert np.allclose(loglik, expected)