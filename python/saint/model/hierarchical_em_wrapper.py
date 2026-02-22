# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 19:10:58 2026

@author: Erik
"""

# %% Imports

import numpy as np


# %% Hierarchical EM wrapper

def run_em_hierarchical(
    X,
    hyperparams,
    biological_bait,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False
):
    """
    Run the hierarchical EM algorithm for a single bait. This function takes an n by r
    matrix of replicate level spectral counts and performs alternating E and M steps
    to update mixture component parameters for each prey under a hierarchical model
    with Dirichlet and Gamma priors estimated by empirical Bayes.

    X is an n by r matrix of replicate level spectral counts.
    hyperparams contains initial values for lambda1, lambda2, lambda3, and pi.
    biological_bait is the biological name associated with this bait.
    max_iter, tol_loglik, tol_params, seed, and verbose control EM behavior.

    The output is a dictionary containing histories of loglik, lambda1, lambda2,
    lambda3, tau, pi, and gamma, as well as their final values and prior histories.
    """

    np.random.seed(seed)

    n, r = X.shape

    lambda1 = hyperparams["lambda1_init"].astype(float).copy()
    lambda2 = hyperparams["lambda2_init"].astype(float).copy()
    lambda3 = hyperparams["lambda3_init"].astype(float).copy()

    tau = np.ones(n)

    pi = hyperparams["pi_init"].astype(float).copy()

    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    lambda3_history = []
    tau_history = []
    pi_history = []
    gamma_history = []

    alpha_history = []
    a_history = []
    b_history = []

    X_sum = X.sum(axis=1)
    eps = 1e-12

    alpha = np.ones(3, dtype=float) * 2.0

    def gamma_moments_from_lambda(lam):
        m = float(np.mean(lam))
        v = float(np.var(lam))
        v = max(v, eps)
        a = m * m / v
        b = m / v
        a = max(a, 1.0)
        b = max(b, eps)
        return a, b

    a1, b1 = gamma_moments_from_lambda(lambda1)
    a2, b2 = gamma_moments_from_lambda(lambda2)
    a3, b3 = gamma_moments_from_lambda(lambda3)

    for it in range(max_iter):

        # %% E step

        log_lambda1 = np.log(lambda1 + eps)
        log_lambda2 = np.log(lambda2 + eps)
        log_lambda3 = np.log(lambda3 + eps)

        loglik1 = X_sum * log_lambda1 - r * lambda1
        loglik2 = X_sum * log_lambda2 - r * lambda2
        loglik3 = X_sum * log_lambda3 - r * lambda3

        log_pi = np.log(pi + eps)

        log_post1 = log_pi[0] + loglik1
        log_post2 = log_pi[1] + loglik2
        log_post3 = log_pi[2] + loglik3

        log_posts = np.vstack([log_post1, log_post2, log_post3]).T
        max_log = np.max(log_posts, axis=1, keepdims=True)
        stabilized = log_posts - max_log
        exp_posts = np.exp(stabilized)
        denom = exp_posts.sum(axis=1, keepdims=True)
        gamma = exp_posts / denom

        gamma1 = gamma[:, 0]
        gamma2 = gamma[:, 1]
        gamma3 = gamma[:, 2]

        loglik = float(np.sum(max_log[:, 0] + np.log(denom[:, 0] + eps)))

        # %% Empirical Bayes updates for priors

        m = gamma.mean(axis=0)
        v = gamma.var(axis=0) + eps
        alpha0 = (m * (1.0 - m) / v - 1.0).clip(min=1.0)
        alpha = m * alpha0
        alpha = np.maximum(alpha, 1.0)

        a1, b1 = gamma_moments_from_lambda(lambda1)
        a2, b2 = gamma_moments_from_lambda(lambda2)
        a3, b3 = gamma_moments_from_lambda(lambda3)

        # %% M step

        num1 = gamma1 * X_sum + (a1 - 1.0)
        num2 = gamma2 * X_sum + (a2 - 1.0)
        num3 = gamma3 * X_sum + (a3 - 1.0)

        den1 = gamma1 * r + b1 + eps
        den2 = gamma2 * r + b2 + eps
        den3 = gamma3 * r + b3 + eps

        lambda1_new = num1 / den1
        lambda2_new = num2 / den2
        lambda3_new = num3 / den3

        lambda1_new = np.clip(lambda1_new, 1e-6, 1e6)
        lambda2_new = np.clip(lambda2_new, 1e-6, 1e6)
        lambda3_new = np.clip(lambda3_new, 1e-6, 1e6)

        pi_new = (gamma.sum(axis=0) + alpha - 1.0)
        pi_new = pi_new / np.sum(pi_new)

        loglik_history.append(loglik)
        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        lambda3_history.append(lambda3_new.copy())
        tau_history.append(tau.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())

        alpha_history.append(alpha.copy())
        a_history.append(np.array([a1, a2, a3], dtype=float))
        b_history.append(np.array([b1, b2, b3], dtype=float))

        delta_lambda = (
            np.max(np.abs(lambda1_new - lambda1)) +
            np.max(np.abs(lambda2_new - lambda2)) +
            np.max(np.abs(lambda3_new - lambda3))
        )
        delta_pi = np.max(np.abs(pi_new - pi))

        if verbose:
            print(
                f"Iter {it} loglik {loglik:.4f} "
                f"delta_lambda {delta_lambda:.4e} delta_pi {delta_pi:.4e}"
            )

        lambda1 = lambda1_new
        lambda2 = lambda2_new
        lambda3 = lambda3_new
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
        "lambda3_history": lambda3_history,
        "tau_history": tau_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history,
        "alpha_history": alpha_history,
        "a_history": a_history,
        "b_history": b_history,
        "lambda1": lambda1,
        "lambda2": lambda2,
        "lambda3": lambda3,
        "tau": tau,
        "pi": pi,
        "gamma": gamma
    }

