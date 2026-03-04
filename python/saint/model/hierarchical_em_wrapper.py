# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:08:34 2026

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
    verbose=False,
):
    """
    Hierarchical EM algorithm for a single bait under a three-component Poisson
    mixture model with empirical Bayes updates for the prior hyperparameters.

    Model structure

    For a given bait, X is an n × r matrix of spectral counts:
      n: number of prey proteins
      r: number of technical or biological replicates

    Each prey belongs to one of three latent components:

      • Component 1 (index 0): background or non-interactor
        Poisson rate: lambda1

      • Component 2 (index 1): contaminant or high background
        Poisson rate: lambda2

      • Component 3 (index 2): true signal or interactor
        Poisson rate: lambda3

    Priors

    Mixture weights pi follow a Dirichlet distribution with parameters alpha.
    Each lambda_k follows a Gamma(a_k, b_k) prior. Hyperparameters are updated
    via empirical Bayes moment-matching based on the current lambda estimates.

    EM algorithm

    E step:
        Compute gamma[i, k], the posterior probability that prey i belongs to
        component k.

    M step:
        Update lambda1, lambda2, lambda3, and pi by maximizing the expected
        posterior objective, combining Poisson likelihood terms weighted by
        gamma and the contributions from the Dirichlet and Gamma priors.

    Returns

    A dictionary containing:
        loglik_history
        lambda1_history, lambda2_history, lambda3_history
        tau_history
        pi_history
        gamma_history
        alpha_history
        a_history, b_history
        lambda1, lambda2, lambda3, tau, pi, gamma

    Parameters

    X : array
        Data matrix for the current bait with shape n × r.
    hyperparams : dict
        Optional hyperparameters for the hierarchical EM model.
    biological_bait : str
        Identifier for the current bait.
    max_iter : int
        Maximum number of EM iterations.
    tol_loglik : float
        Convergence tolerance for changes in the log likelihood.
    tol_params : float
        Convergence tolerance for changes in parameter estimates.
    seed : int
        Random seed used for initialization.
    verbose : bool
        If True, print progress information during EM.
    """

    # %% Initialization and set-up
    if hyperparams is None:
        hyperparams = {}

    # Dirichlet prior defaults
    hyperparams.setdefault("alpha1", 2.0)
    hyperparams.setdefault("alpha2", 2.0)
    hyperparams.setdefault("alpha3", 2.0)

    # Mixing proportions initialization (uniform)
    hyperparams.setdefault("pi_init", np.ones(3, dtype=float) / 3.0)

    # Global shrinkage parameter
    hyperparams.setdefault("tau", 1.0)

    # Tau grid (not used here but kept for API stability)
    hyperparams.setdefault("tau_grid", [0.1, 0.5, 1.0, 2.0])

    np.random.seed(seed)

    # Dimensions
    n, r = X.shape

    # Prey-dimensional lambda initializations
    global_mean = X.mean() / r

    lambda1 = np.full(n, 0.8 * global_mean)
    lambda2 = np.full(n, 1.0 * global_mean)
    lambda3 = np.full(n, 1.2 * global_mean)


    # Global shrinkage strength
    tau = float(hyperparams.get("tau", 1.0))

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
    alpha = np.array(
        [hyperparams["alpha1"], hyperparams["alpha2"], hyperparams["alpha3"]],
        dtype=float,
    )

    # Empirical Bayes hyperparameter update
    def gamma_moments_from_lambda(lam):
        m = float(np.mean(lam))
        v = float(np.var(lam))
        v = max(v, eps)
        a = max(m * m / v, 1.0)
        b = max(m / v, eps)
        return a, b

    # Initial Gamma hyperparameters
    a1, b1 = gamma_moments_from_lambda(lambda1)
    a2, b2 = gamma_moments_from_lambda(lambda2)
    a3, b3 = gamma_moments_from_lambda(lambda3)
    
    # %% EM iterations
    for it in range(max_iter):
    
        # %% E-step
    
        # Component-wise log lambda
        log_lambda1 = np.log(lambda1 + eps)
        log_lambda2 = np.log(lambda2 + eps)
        log_lambda3 = np.log(lambda3 + eps)
    
        # Component-wise log-likelihood contributions: log p(X_i | lambda_k)
        loglik1 = X_sum * log_lambda1 - r * lambda1
        loglik2 = X_sum * log_lambda2 - r * lambda2
        loglik3 = X_sum * log_lambda3 - r * lambda3
    
        # Log mixture weights
        log_pi = np.log(pi + eps)
    
        # Posterior log-likelihood for each component
        log_post1 = log_pi[0] + loglik1
        log_post2 = log_pi[1] + loglik2
        log_post3 = log_pi[2] + loglik3
    
        # Shape (n, 3)
        posterior_loglik = np.vstack([log_post1, log_post2, log_post3]).T
    
        # Stabilize before exponentiating
        max_log = np.max(posterior_loglik, axis=1, keepdims=True)
        stabilized_log_posterior = posterior_loglik - max_log
    
        # Convert posterior log-likelihood to posterior weights
        posterior_weights = np.exp(stabilized_log_posterior)
    
        # Normalizing constant
        normalizing_constant = posterior_weights.sum(axis=1, keepdims=True)
    
        # Posterior responsibilities
        gamma = posterior_weights / (normalizing_constant + eps)
    
        gamma1 = gamma[:, 0]
        gamma2 = gamma[:, 1]
        gamma3 = gamma[:, 2]
    
        # Log-likelihood for this iteration
        loglik = float(
            np.sum(max_log[:, 0] + np.log(normalizing_constant[:, 0] + eps))
        )
    
        # %% Empirical Bayes updates
    
        # Update Dirichlet prior on mixture weights
        m = gamma.mean(axis=0)
        v = gamma.var(axis=0) + eps
        alpha0 = (m * (1.0 - m) / v - 1.0).clip(min=1.0)
        alpha = np.maximum(m * alpha0, 1.0)
    
        # Gamma(a, b) hyperparameters for each component
        # (moment-matching on lambda)
        # NOTE: These are intentionally *not* updated each iteration
        # unless explicitly desired. The original code commented them out.
        # We preserve that behavior for stability.
        # a1, b1 = gamma_moments_from_lambda(lambda1)
        # a2, b2 = gamma_moments_from_lambda(lambda2)
        # a3, b3 = gamma_moments_from_lambda(lambda3)
    
        # %% M-step
    
        # Component 1
        data_contrib_1 = gamma1 * X_sum
        prior_contrib_1 = tau * (a1 - 1.0)
        data_weight_1 = gamma1 * r
        prior_weight_1 = tau * b1
    
        lambda1_new = np.clip(
            (data_contrib_1 + prior_contrib_1) /
            (data_weight_1 + prior_weight_1 + eps),
            1e-6, 1e6
        )
    
        # Component 2
        data_contrib_2 = gamma2 * X_sum
        prior_contrib_2 = tau * (a2 - 1.0)
        data_weight_2 = gamma2 * r
        prior_weight_2 = tau * b2
    
        lambda2_new = np.clip(
            (data_contrib_2 + prior_contrib_2) /
            (data_weight_2 + prior_weight_2 + eps),
            1e-6, 1e6
        )
    
        # Component 3
        data_contrib_3 = gamma3 * X_sum
        prior_contrib_3 = tau * (a3 - 1.0)
        data_weight_3 = gamma3 * r
        prior_weight_3 = tau * b3
    
        lambda3_new = np.clip(
            (data_contrib_3 + prior_contrib_3) /
            (data_weight_3 + prior_weight_3 + eps),
            1e-6, 1e6
        )
    
        # Updated mixing proportions
        pi_new = gamma.sum(axis=0) + alpha - 1.0
        pi_new = pi_new / np.sum(pi_new)
    
        # Store histories
        loglik_history.append(loglik)
        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        lambda3_history.append(lambda3_new.copy())
        tau_history.append(tau)
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())
    
        alpha_history.append(alpha.copy())
        a_history.append(np.array([a1, a2, a3], dtype=float))
        b_history.append(np.array([b1, b2, b3], dtype=float))
    
        # Convergence diagnostics
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
        "lambda1": lambda1.copy(),
        "lambda2": lambda2.copy(),
        "lambda3": lambda3.copy(),
        "tau": float(tau),
        "pi": pi.copy(),
        "gamma": gamma.copy(),
    }