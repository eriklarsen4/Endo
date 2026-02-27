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
    Hierarchical EM algorithm for a single bait under a three‑component mixture model
    with empirical Bayes updates for the priors.

    Model structure
    ---------------
    For a given bait, we observe an n × r matrix X of spectral counts, where:
      - n: number of prey proteins
      - r: number of technical/biological replicates for this bait

    Each prey is modeled as arising from one of three latent components:

      • Component 1 (index 0): background / non‑interactor
        - Poisson rate parameter: lambda1
        - Represents proteins with only stochastic or technical counts.

      • Component 2 (index 1): true signal / interactor
        - Poisson rate parameter: lambda2
        - Represents proteins enriched due to true biological interaction.

      • Component 3 (index 2): contaminant / high‑background
        - Poisson rate parameter: lambda3
        - Represents proteins with elevated counts not attributable to specific
          interaction (e.g., sticky binders, common contaminants).

    The model is hierarchical in two ways:
      1. Mixture weights pi have a Dirichlet prior with parameter alpha.
      2. Each lambda_k (k = 1, 2, 3) has a Gamma(a_k, b_k) prior.
         The hyperparameters (alpha, a_k, b_k) are updated via empirical Bayes
         using moment‑matching based on the current parameter estimates.

    EM algorithm overview
    ---------------------
    The EM algorithm alternates between:

      • E‑step:
          Computes gamma[i, k], the posterior probability that prey i belongs
          to mixture component k, integrating over replicate counts via X_sum.
          These probabilities represent the expected values of the unobserved
          latent class indicators Z_{ik}.

      • M‑step:
          Updates lambda1, lambda2, lambda3, and pi by maximizing the expected
          log‑likelihood of the *complete‑data hierarchical model*. Here,
          “complete data” refers specifically to:
              (a) the observed replicate‑level counts X (or their row sums X_sum),
              (b) the unobserved latent component assignments Z, and
              (c) the latent rates lambda_k under their Gamma priors.
          The M‑step maximizes:
              E_{Z, lambda}[ log p(X, Z, lambda, pi | alpha, a, b) ],
          where the expectation is taken with respect to the posterior
          distributions of Z and lambda implied by the current parameters.

    Parameter definitions
    ---------------------
    X : array of shape (n, r)
        Replicate‑level spectral counts for n prey proteins and r replicates.

    hyperparams : dict
        Contains initial values:
            "lambda1_init" : initial Poisson rate for background component
            "lambda2_init" : initial Poisson rate for signal component
            "lambda3_init" : initial Poisson rate for contaminant component
            "pi_init"      : initial mixing proportions (length 3)

    biological_bait : str
        Biological name of the bait (used for labeling downstream results).

    max_iter : int
        Maximum number of EM iterations.

    tol_loglik : float
        Convergence tolerance for log‑likelihood change.

    tol_params : float
        Convergence tolerance for parameter changes.

    seed : int
        Random seed for reproducibility.

    verbose : bool
        If True, prints iteration‑level diagnostics.

    Returns
    -------
    dict containing:
        loglik_history : list of float
        lambda1_history, lambda2_history, lambda3_history : list of arrays
        tau_history : list of arrays (currently placeholder, kept for API stability)
        pi_history : list of arrays
        gamma_history : list of (n × 3) arrays
        alpha_history : list of Dirichlet parameter vectors
        a_history, b_history : list of Gamma hyperparameter vectors
        lambda1, lambda2, lambda3, tau, pi, gamma : final parameter values

    Notes
    -----
    Variable names (lambda1, lambda2, lambda3, pi, gamma, alpha, a, b, tau)
    are kept for architectural compatibility with downstream code, but each
    is documented explicitly to clarify its biological and statistical meaning.
    """
    
    # %% Initialization and set-up
    # Set random seed for reproducibility of any randomized initialization
    np.random.seed(seed)
    
    # Dimensions of the replicate-level data matrix
    # n: number of prey proteins
    # r: number of replicates for this bait
    n, r = X.shape
    
    # Initialize Poisson rate parameters for the three mixture components.
    # lambda1: background / non-interactor rate (vector of length n)
    # lambda2: noise / contaminant rate (vector of length n)
    # lambda3: true signal / interactor rate (vector of length n)
    lambda1 = hyperparams["lambda1_init"].astype(float).copy()
    lambda2 = hyperparams["lambda2_init"].astype(float).copy()
    lambda3 = hyperparams["lambda3_init"].astype(float).copy()
    
    # tau is currently a placeholder for potential prey-level random effects.
    # It is kept for API stability and future extensibility.
    tau = np.ones(n)
    
    # Initialize mixing proportions for the three components.
    # pi[0]: background
    # pi[1]: noise / contaminant
    # pi[2]: true signal
    pi = hyperparams["pi_init"].astype(float).copy()
    
    # Histories for monitoring EM convergence and parameter trajectories
    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    lambda3_history = []
    tau_history = []
    pi_history = []
    gamma_history = []
    
    # Histories for empirical Bayes hyperparameters
    alpha_history = []   # Dirichlet parameters for pi
    a_history = []       # Gamma shape parameters for lambda_k
    b_history = []       # Gamma rate parameters for lambda_k
    
    # Precompute row sums across replicates for Poisson likelihood
    X_sum = X.sum(axis=1)
    
    # Small constant to avoid log(0) and division by zero
    eps = 1e-12
    
    # Initialize Dirichlet prior parameters for pi
    # alpha[k] >= 1 ensures a proper, non-degenerate prior
    alpha = np.ones(3, dtype=float) * 2.0
    
    # Helper: moment-matching to estimate Gamma(a, b) hyperparameters
    # for each lambda_k vector.
    def gamma_moments_from_lambda(lam):
        # Mean and variance of current lambda estimates
        m = float(np.mean(lam))
        v = float(np.var(lam))
        v = max(v, eps)
    
        # Moment-matching for Gamma(a, b):
        #   mean = a / b  =>  a = m^2 / v
        #   var  = a / b^2 => b = m / v
        a = m * m / v
        b = m / v
    
        # Enforce minimal shape and rate to avoid degeneracy
        a = max(a, 1.0)
        b = max(b, eps)
        return a, b
    
    # Initial Gamma hyperparameters for each component
    a1, b1 = gamma_moments_from_lambda(lambda1)
    a2, b2 = gamma_moments_from_lambda(lambda2)
    a3, b3 = gamma_moments_from_lambda(lambda3)
    
    
    
    for it in range(max_iter):

    # %% E-step: posterior membership probabilities gamma[i, k]
    # for each prey i and component k (background, contaminant, signal).

    # Log of Poisson rates for numerical stability
    log_lambda1 = np.log(lambda1 + eps)
    log_lambda2 = np.log(lambda2 + eps)
    log_lambda3 = np.log(lambda3 + eps)

    # Component-specific log-likelihood contributions:
    # log p(X_sum[i] | lambda_k) up to an additive constant.
    # Note: replicate count r multiplies the rate term.
    loglik1 = X_sum * log_lambda1 - r * lambda1
    loglik2 = X_sum * log_lambda2 - r * lambda2
    loglik3 = X_sum * log_lambda3 - r * lambda3

    # Log mixing proportions
    log_pi = np.log(pi + eps)

    # Unnormalized log posterior for each component
    log_post1 = log_pi[0] + loglik1
    log_post2 = log_pi[1] + loglik2
    log_post3 = log_pi[2] + loglik3

    # Stack into an (n × 3) matrix of log posteriors
    log_posts = np.vstack([log_post1, log_post2, log_post3]).T

    # Stabilize with log-sum-exp trick
    max_log = np.max(log_posts, axis=1, keepdims=True)
    stabilized = log_posts - max_log
    exp_posts = np.exp(stabilized)
    denom = exp_posts.sum(axis=1, keepdims=True)

    # Posterior membership probabilities gamma[i, k]
    gamma = exp_posts / denom

    gamma1 = gamma[:, 0]   # background responsibility
    gamma2 = gamma[:, 1]   # contaminant responsibility
    gamma3 = gamma[:, 2]   # signal responsibility

    # Observed-data log-likelihood
    loglik = float(np.sum(max_log[:, 0] + np.log(denom[:, 0] + eps)))

    # %% Empirical Bayes updates for priors
    # Update Dirichlet prior alpha for pi using moment-matching on gamma.
    m = gamma.mean(axis=0)
    v = gamma.var(axis=0) + eps
    alpha0 = (m * (1.0 - m) / v - 1.0).clip(min=1.0)
    alpha = m * alpha0
    alpha = np.maximum(alpha, 1.0)

    # Update Gamma(a_k, b_k) hyperparameters for each lambda_k
    a1, b1 = gamma_moments_from_lambda(lambda1)
    a2, b2 = gamma_moments_from_lambda(lambda2)
    a3, b3 = gamma_moments_from_lambda(lambda3)

    # %% M-step: update lambda1, lambda2, lambda3, and pi

    # Weighted sufficient statistics for Poisson rates with Gamma prior
    num1 = gamma1 * X_sum + (a1 - 1.0)
    num2 = gamma2 * X_sum + (a2 - 1.0)
    num3 = gamma3 * X_sum + (a3 - 1.0)

    den1 = gamma1 * r + b1 + eps
    den2 = gamma2 * r + b2 + eps
    den3 = gamma3 * r + b3 + eps

    lambda1_new = num1 / den1
    lambda2_new = num2 / den2
    lambda3_new = num3 / den3

    # Numerical safety
    lambda1_new = np.clip(lambda1_new, 1e-6, 1e6)
    lambda2_new = np.clip(lambda2_new, 1e-6, 1e6)
    lambda3_new = np.clip(lambda3_new, 1e-6, 1e6)

    # Updated mixing proportions with Dirichlet prior
    pi_new = (gamma.sum(axis=0) + alpha - 1.0)
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

    # Convergence diagnostics
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

    # Commit updates
    lambda1 = lambda1_new
    lambda2 = lambda2_new
    lambda3 = lambda3_new
    pi = pi_new

    # %% Convergence check
    if it > 0:
        if (
            abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik
            and delta_lambda < tol_params
            and delta_pi < tol_params
        ):
            break


    # %% Final output: histories and last-iteration parameter values
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