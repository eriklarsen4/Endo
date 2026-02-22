# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 19:09:57 2026

@author: Erik
"""

# %% Imports

import numpy as np


# %% Classical EM wrapper

def run_em_classical(
    X_sum,
    hyperparams,
    biological_bait,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False
):
    """
    Run the classical EM algorithm for a single bait. This function takes a vector of
    collapsed spectral counts and performs alternating E and M steps to update mixture
    component parameters for each prey under the classical two component model.

    X_sum is a length n vector of collapsed spectral counts.
    hyperparams contains initial values for lambda1, lambda2, and pi.
    biological_bait is the biological name associated with this bait.
    max_iter, tol_loglik, tol_params, seed, and verbose control EM behavior.

    The output is a dictionary containing histories of loglik, lambda1, lambda2, pi,
    gamma, as well as their final values.
    """

    np.random.seed(seed)

    n = X_sum.shape[0]

    lambda1 = hyperparams["lambda1_init"].astype(float).copy()
    lambda2 = hyperparams["lambda2_init"].astype(float).copy()
    pi = hyperparams["pi_init"].astype(float).copy()

    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    pi_history = []
    gamma_history = []

    eps = 1e-12

    for it in range(max_iter):

        # %% E step

        log_lambda1 = np.log(lambda1 + eps)
        log_lambda2 = np.log(lambda2 + eps)

        loglik1 = X_sum * log_lambda1 - lambda1
        loglik2 = X_sum * log_lambda2 - lambda2

        log_pi = np.log(pi + eps)

        log_post1 = log_pi[0] + loglik1
        log_post2 = log_pi[1] + loglik2

        log_posts = np.vstack([log_post1, log_post2]).T
        max_log = np.max(log_posts, axis=1, keepdims=True)
        stabilized = log_posts - max_log
        exp_posts = np.exp(stabilized)
        denom = exp_posts.sum(axis=1, keepdims=True)
        gamma = exp_posts / denom

        gamma1 = gamma[:, 0]
        gamma2 = gamma[:, 1]

        # log likelihood using log sum exp
        loglik = float(np.sum(max_log[:, 0] + np.log(denom[:, 0] + eps)))

        # %% M step

        num1 = gamma1 * X_sum
        num2 = gamma2 * X_sum

        den1 = gamma1 + eps
        den2 = gamma2 + eps

        lambda1_new = num1 / den1
        lambda2_new = num2 / den2

        lambda1_new = np.clip(lambda1_new, 1e-6, 1e6)
        lambda2_new = np.clip(lambda2_new, 1e-6, 1e6)

        pi_new = np.array([
            gamma1.sum() / n,
            gamma2.sum() / n
        ], dtype=float)

        loglik_history.append(loglik)
        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())

        delta_lambda = (
            np.max(np.abs(lambda1_new - lambda1)) +
            np.max(np.abs(lambda2_new - lambda2))
        )
        delta_pi = np.max(np.abs(pi_new - pi))

        if verbose:
            print(
                f"Iter {it} loglik {loglik:.4f} "
                f"delta_lambda {delta_lambda:.4e} delta_pi {delta_pi:.4e}"
            )

        lambda1 = lambda1_new
        lambda2 = lambda2_new
        pi = pi_new

        if it > 0:
            if (
                abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik
                and delta_lambda < tol_params
                and delta_pi < tol_params
            ):
                break

    return {
        "loglik_history": loglik_history,
        "lambda1_history": lambda1_history,
        "lambda2_history": lambda2_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history,
        "lambda1": lambda1,
        "lambda2": lambda2,
        "pi": pi,
        "gamma": gamma
    }

