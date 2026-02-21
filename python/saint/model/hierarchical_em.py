# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 08:18:37 2026

@author: Erik
"""

# %% Hierarchical SAINT EM

import numpy as np

from .responsibilities import compute_responsibilities
from .updates.lambda_updates import update_lambda_hierarchical
from .updates.pi_updates import update_pi
from .init_tau import init_tau
from .likelihood import compute_loglik


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
    lambda3 = hyperparams["lambda3_init"].copy()
    pi = hyperparams["pi_init"].copy()

    if "alpha" in hyperparams:
        alpha = hyperparams["alpha"].copy()
    else:
        alpha = np.ones(3, dtype=float)

    tau = init_tau()

    loglik_hist = []
    lambda1_hist = []
    lambda2_hist = []
    lambda3_hist = []
    pi_hist = []
    gamma_hist = []
    tau_hist = []

    prev_loglik = None

    for iteration in range(max_iter):

        gamma = compute_responsibilities(
            X,
            [lambda1, lambda2, lambda3],
            pi
        )

        lambda1_new, lambda2_new, lambda3_new = update_lambda_hierarchical(
            X,
            gamma,
            alpha,
            tau
        )

        pi_new = update_pi(gamma)

        tau_new = tau  # placeholder if tau is later updated

        loglik = compute_loglik(
            X,
            [lambda1_new, lambda2_new, lambda3_new],
            pi_new
        )

        loglik_hist.append(loglik)
        lambda1_hist.append(lambda1_new.copy())
        lambda2_hist.append(lambda2_new.copy())
        lambda3_hist.append(lambda3_new.copy())
        pi_hist.append(pi_new.copy())
        gamma_hist.append(gamma.copy())
        tau_hist.append(tau_new)

        if verbose:
            print(f"Iteration {iteration}, loglik = {loglik:.6f}")

        if prev_loglik is not None:
            if abs(loglik - prev_loglik) < tol_loglik:
                break

            delta = (
                np.sum(np.abs(lambda1_new - lambda1))
                + np.sum(np.abs(lambda2_new - lambda2))
                + np.sum(np.abs(lambda3_new - lambda3))
                + np.sum(np.abs(pi_new - pi))
            )
            if delta < tol_params:
                break

        lambda1 = lambda1_new
        lambda2 = lambda2_new
        lambda3 = lambda3_new
        pi = pi_new
        tau = tau_new
        prev_loglik = loglik

    return {
        "loglik_history": loglik_hist,
        "lambda1_history": lambda1_hist,
        "lambda2_history": lambda2_hist,
        "lambda3_history": lambda3_hist,
        "pi_history": pi_hist,
        "gamma_history": gamma_hist,
        "tau_history": tau_hist
    }



