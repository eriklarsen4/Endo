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
    Hierarchical EM algorithm for a single bait under a three‑component mixture
    model with empirical Bayes updates for the priors.

    Model structure
    ---------------
    For a given bait, we observe an n × r matrix X of spectral counts:
      - n: number of prey proteins
      - r: number of technical/biological replicates

    Each prey belongs to one of three latent components:

      • Component 1 (index 0): background / non‑interactor
        - Poisson rate: lambda1
        - Represents proteins with only stochastic or technical counts.

      • Component 2 (index 1): contaminant / high‑background
        - Poisson rate: lambda2
        - Represents proteins with elevated counts not attributable to specific
          interaction (sticky binders, common contaminants).

      • Component 3 (index 2): true signal / interactor
        - Poisson rate: lambda3
        - Represents proteins enriched due to true biological interaction.

    Priors:
      - Mixture weights pi ~ Dirichlet(alpha)
      - Each lambda_k ~ Gamma(a_k, b_k)
      - Hyperparameters (alpha, a_k, b_k) updated via empirical Bayes.

    EM algorithm
    ------------
    E‑step:
        Compute gamma[i, k], the posterior probability that prey i belongs to
        component k, using X_sum = row sums of X.

    M‑step:
        Update lambda1, lambda2, lambda3, and pi by maximizing the expected
        complete‑data log‑likelihood under the hierarchical model.

    Returns
    -------
    A dictionary containing:
        loglik_history
        lambda1_history, lambda2_history, lambda3_history
        tau_history
        pi_history
        gamma_history
        alpha_history
        a_history, b_history
        lambda1, lambda2, lambda3, tau, pi, gamma
    """

    # %% Initialization and set-up
    np.random.seed(seed)

    # Dimensions
    n, r = X.shape

    # ORIGINAL INITIALIZATIONS (preserved exactly)
    lambda1 = X.mean(axis=1).astype(float).copy()
    lambda2 = np.full(n, 0.5, dtype=float)
    lambda3 = np.full(n, 0.1, dtype=float)

    # tau placeholder
    tau = np.ones(n)

    # Mixing proportions
    pi = hyperparams["pi_init"].astype(float).copy()

    # Histories
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

    # Precompute row sums
    X_sum = X.sum(axis=1)

    eps = 1e-12

    # Dirichlet prior for pi
    alpha = np.ones(3, dtype=float) * 2.0

    # Gamma moment-matching helper
    def gamma_moments_from_lambda(lam):
        m = float(np.mean(lam))
        v = float(np.var(lam))
        v = max(v, eps)
        a = m * m / v
        b = m / v
        a = max(a, 1.0)
        b = max(b, eps)
        return a, b

    # Initial Gamma hyperparameters
    a1, b1 = gamma_moments_from_lambda(lambda1)
    a2, b2 = gamma_moments_from_lambda(lambda2)
    a3, b3 = gamma_moments_from_lambda(lambda3)

    # %% EM iterations
    for it in range(max_iter):

        # %% E-step
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

        # %% Empirical Bayes updates
        m = gamma.mean(axis=0)
        v = gamma.var(axis=0) + eps
        alpha0 = (m * (1.0 - m) / v - 1.0).clip(min=1.0)
        alpha = np.maximum(m * alpha0, 1.0)

        a1, b1 = gamma_moments_from_lambda(lambda1)
        a2, b2 = gamma_moments_from_lambda(lambda2)
        a3, b3 = gamma_moments_from_lambda(lambda3)

        # %% M-step
        num1 = gamma1 * X_sum + (a1 - 1.0)
        num2 = gamma2 * X_sum + (a2 - 1.0)
        num3 = gamma3 * X_sum + (a3 - 1.0)

        den1 = gamma1 * r + b1 + eps
        den2 = gamma2 * r + b2 + eps
        den3 = gamma3 * r + b3 + eps

        lambda1_new = np.clip(num1 / den1, 1e-6, 1e6)
        lambda2_new = np.clip(num2 / den2, 1e-6, 1e6)
        lambda3_new = np.clip(num3 / den3, 1e-6, 1e6)

        pi_new = gamma.sum(axis=0) + alpha - 1.0
        pi_new = pi_new / np.sum(pi_new)

        # Store histories
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

        # Convergence
        delta_lambda = (
            np.max(np.abs(lambda1_new - lambda1)) +
            np.max(np.abs(lambda2_new - lambda2)) +
            np.max(np.abs(lambda3_new - lambda3))
        )
        delta_pi = np.max(np.abs(pi_new - pi))

        if it > 0:
            if (
                abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik
                and delta_lambda < tol_params
                and delta_pi < tol_params
            ):
                break

        # Commit updates
        lambda1 = lambda1_new
        lambda2 = lambda2_new
        lambda3 = lambda3_new
        pi = pi_new

    # %% Final output
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