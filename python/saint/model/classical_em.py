# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 00:05:05 2026

@author: Erik
"""

# %% Imports

import numpy as np

from saint.model.responsibilities import compute_responsibilities
from saint.model.updates.lambda_updates import update_lambda_placeholder
from saint.model.updates.pi_updates import update_pi
from saint.model.loglik import compute_loglik


# %% Classical SAINT EM

def run_em(
    X,
    hyperparams,
    max_iter=100,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=None,
    verbose=False
):
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

        gamma = compute_responsibilities(X, lambda1, lambda2, pi)

        lambda_updates = update_lambda_placeholder(X, gamma)
        lambda1_new = lambda_updates[0]
        lambda2_new = lambda_updates[1]

        pi_new = update_pi(gamma)

        loglik = compute_loglik(X, lambda1_new, lambda2_new, pi_new)
        loglik_history.append(loglik)

        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())

        if verbose:
            print(f"Iteration {iteration}, loglik = {loglik:.6f}")

        if prev_loglik is not None:
            if abs(loglik - prev_loglik) < tol_loglik:
                break

            delta = (
                np.sum(np.abs(lambda1_new - lambda1))
                + np.sum(np.abs(lambda2_new - lambda2))
                + np.sum(np.abs(pi_new - pi))
            )
            if delta < tol_params:
                break

        lambda1 = lambda1_new
        lambda2 = lambda2_new
        pi = pi_new
        prev_loglik = loglik

    return {
        "loglik_history": loglik_history,
        "lambda1_history": lambda1_history,
        "lambda2_history": lambda2_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history
    }

