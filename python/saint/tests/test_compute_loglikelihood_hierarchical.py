# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 09:14:57 2026

@author: Erik
"""

# %% Imports
import numpy as np
from scipy.special import gammaln

from saint.model.hierarchical_likelihood import compute_loglik


# %%
def test_compute_hierarchical_loglik_simple_case():
    """
    Minimal test verifying that compute_hierarchical_loglik matches a manually
    computed mixture-of-Poissons likelihood on a tiny dataset, ignoring priors.
    """
    X = np.array([
        [1, 2],
        [0, 1],
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    lambda3 = np.array([4.0, 4.0])

    pi = np.array([0.3, 0.3, 0.4])

    # manual expected value
    expected = 0.0
    for p in range(X.shape[0]):
        Xp = X[p]

        ll1 = (Xp * np.log(lambda1[p]) - lambda1[p] - gammaln(Xp + 1)).sum()
        ll2 = (Xp * np.log(lambda2[p]) - lambda2[p] - gammaln(Xp + 1)).sum()
        ll3 = (Xp * np.log(lambda3[p]) - lambda3[p] - gammaln(Xp + 1)).sum()

        joint = (
            pi[0] * np.exp(ll1) +
            pi[1] * np.exp(ll2) +
            pi[2] * np.exp(ll3)
        )
        expected += np.log(joint)

    loglik = compute_loglik(
        X,
        lambda1,
        lambda2,
        lambda3,
        pi,
        use_exact_poisson_likelihood=True,
        include_priors=False,  # or whatever your signature uses
    )

    assert np.isclose(loglik, expected)

