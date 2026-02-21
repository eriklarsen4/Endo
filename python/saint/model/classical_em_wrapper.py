# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:50:39 2026

@author: Erik
"""

# %% Imports

import numpy as np
from saint.model.classical_responsibilities import compute_responsibilities
from saint.model.updates.lambda_updates import update_lambda1, update_lambda2
from saint.model.updates.pi_updates import update_pi
from saint.model.likelihood import compute_loglik


# %% EM wrapper

def run_em(
    X,
    hyperparams,
    biological_bait,
    max_iter=100,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=None,
    verbose=False
):
    """
    Classical EM algorithm for SAINT.

    lambda1 is the background Poisson rate for each prey (protein).
    lambda2 is the signal Poisson rate for each prey (protein).
    pi is the vector of mixture proportions for the two components.
    gamma represents the posterior probability that each prey belongs
    to each component.

    The E step computes gamma using the current lambda and pi values.
    The M step updates lambda1, lambda2, and pi using gamma.

    X is a numeric matrix of prey (protein) by bait.
    hyperparams is a dictionary containing initial lambda and pi values.
    biological_bait is the protein symbol corresponding to the target bait.
    """

    if seed is not None:
        np.random.seed(seed)

    lambda1 = hyperparams["lambda1_init"].copy()
    lambda2 = hyperparams["lambda2_init"].copy()
    pi = hyperparams["pi_init"].copy()

    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    pi_history = []
    gamma_history = []

    prev_loglik = None

    for iteration in range(max_iter):

        # %% E step
        gamma = compute_responsibilities(X, lambda1, lambda2, pi)

        # %% M step
        lambda1_new = update_lambda1(X, gamma)
        lambda2_new = update_lambda2(X, gamma)
        pi_new = update_pi(gamma)

        # %% Convergence check
        loglik = compute_loglik(X, lambda1_new, lambda2_new, pi_new)
        loglik_history.append(loglik)

        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())

        if prev_loglik is not None:
            if abs(loglik - prev_loglik) < tol_loglik:
                break

            if (
                np.max(np.abs(lambda1_new - lambda1)) < tol_params
                and np.max(np.abs(lambda2_new - lambda2)) < tol_params
                and np.max(np.abs(pi_new - pi)) < tol_params
            ):
                break

        lambda1 = lambda1_new
        lambda2 = lambda2_new
        pi = pi_new
        prev_loglik = loglik

        if verbose:
            print("Iteration", iteration, "loglik", loglik)

    return {
        "lambda1_history": lambda1_history,
        "lambda2_history": lambda2_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history,
        "loglik_history": loglik_history,
        "biological_bait": biological_bait
    }

