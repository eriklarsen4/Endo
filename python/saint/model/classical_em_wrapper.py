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
    Classical EM algorithm for a single bait under the two‑component mixture model. Loop over baits.

    Model structure
    ---------------
    Each prey is modeled as arising from one of two latent components:

      • Component 1 (index 0): background / non‑interactor
        - Poisson rate parameter: lambda1
        - Represents proteins with only stochastic or technical spectral counts.

      • Component 2 (index 1): true signal / interactor
        - Poisson rate parameter: lambda2
        - Represents proteins enriched due to true biological interaction.

    EM algorithm overview
    ---------------------
    The EM algorithm alternates between:

      • E‑step:
          Computes gamma[i, k], the posterior probability that prey i belongs
          to mixture component k. These probabilities represent the expected
          values of the unobserved latent class indicators Z_{ik}.

      • M‑step:
          Updates lambda1, lambda2, and pi by maximizing the expected
          log‑likelihood of the *complete‑data model*. Here, “complete data”
          refers specifically to:
              (a) the observed spectral counts X_sum, and
              (b) the unobserved latent component assignments Z.
          The M‑step maximizes:
              E_Z[ log p(X_sum, Z | lambda1, lambda2, pi) ],
          where the expectation is taken with respect to the posterior
          distribution of Z computed in the E‑step.

    Parameter definitions
    ---------------------
    X_sum : array of shape (n,)
        Collapsed spectral counts for n prey proteins for this bait.

    hyperparams : dict
        Contains initial values:
            "lambda1_init" : initial Poisson rate for background component
            "lambda2_init" : initial Poisson rate for signal component
            "pi_init"      : initial mixing proportions (length 2)

    biological_bait : str
        Biological name of the bait (used for labeling downstream results).

    max_iter : int
        Maximum number of EM iterations.

    tol_loglik : float
        Convergence tolerance for log‑likelihood change (the lowest limit of a change in the log-likelihood).

    tol_params : float
        Convergence tolerance for parameter changes (the lowest limit of a change in (any) parameter value).

    seed : int
        Random seed for reproducibility.

    verbose : bool
        If True, prints iteration‑level diagnostics.

    Returns
    -------
    dict containing:
        loglik_history : list of float
        lambda1_history, lambda2_history : list of arrays
        pi_history : list of arrays
        gamma_history : list of (n × 2) arrays
        lambda1, lambda2, pi, gamma : final parameter values

    Notes
    -----
    Variable names (lambda1, lambda2, pi, gamma) are kept for architectural
    compatibility with downstream code, but each is documented explicitly
    to clarify its biological and statistical meaning.
    """

    # %% Initialization and set-up
    # Set random seed for reproducibility of any randomized initialization
    np.random.seed(seed)

    # Number of prey proteins for this bait
    n = X_sum.shape[0]

    # Initialize Poisson rate parameters for the two mixture components
    # lambda1: background / non‑interactor rate (per‑prey, same length as X_sum)
    # lambda2: signal / interactor rate (per‑prey, same length as X_sum)
    lambda1 = hyperparams["lambda1_init"].astype(float).copy()
    lambda2 = hyperparams["lambda2_init"].astype(float).copy()

    # Initialize mixing proportions for the two components:
    # pi[0]: prior probability of background component
    # pi[1]: prior probability of signal component
    pi = hyperparams["pi_init"].astype(float).copy()

    # Histories for monitoring EM convergence and parameter trajectories
    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    pi_history = []
    gamma_history = []

    # Small constant to avoid log(0) and division by zero
    eps = 1e-12

    # Main EM loop: alternate E‑step and M‑step until convergence or max_iter
    for it in range(max_iter):

        # %% E‑step: compute posterior membership probabilities gamma[i, k]
        # for each prey i and component k (background vs signal).

        # Log of Poisson rates for numerical stability
        log_lambda1 = np.log(lambda1 + eps)
        log_lambda2 = np.log(lambda2 + eps)

        # Log‑likelihood contributions for each component:
        # log p(X_sum[i] | lambda_k) up to an additive constant
        loglik1 = X_sum * log_lambda1 - lambda1
        loglik2 = X_sum * log_lambda2 - lambda2

        # Log mixing proportions
        log_pi = np.log(pi + eps)

        # Unnormalized log posterior for each component
        log_post1 = log_pi[0] + loglik1
        log_post2 = log_pi[1] + loglik2

        # Stack into an (n × 2) matrix of log posteriors
        log_posts = np.vstack([log_post1, log_post2]).T

        # Stabilize with log‑sum‑exp trick to avoid numerical underflow
        max_log = np.max(log_posts, axis=1, keepdims=True)
        stabilized = log_posts - max_log
        exp_posts = np.exp(stabilized)
        denom = exp_posts.sum(axis=1, keepdims=True)

        # Posterior membership probabilities gamma[i, k]
        gamma = exp_posts / denom

        # Component‑specific responsibilities
        gamma1 = gamma[:, 0]  # responsibility for background component
        gamma2 = gamma[:, 1]  # responsibility for signal component

        # Observed‑data log‑likelihood using the stabilized log‑sum‑exp form
        loglik = float(np.sum(max_log[:, 0] + np.log(denom[:, 0] + eps)))
        
        # %% M‑step: update lambda1, lambda2, and pi using the posterior weights.
        # Each update corresponds to maximizing the expected complete‑data
        # log‑likelihood under the current responsibilities gamma.
    
        # Weighted sufficient statistics for Poisson rates
        num1 = gamma1 * X_sum          # numerator for background rate update
        num2 = gamma2 * X_sum          # numerator for signal rate update
    
        # Effective sample sizes for each component
        den1 = gamma1 + eps
        den2 = gamma2 + eps
    
        # Updated Poisson rates (element‑wise for each prey)
        lambda1_new = num1 / den1
        lambda2_new = num2 / den2
    
        # Numerical safety: keep rates within a reasonable range
        lambda1_new = np.clip(lambda1_new, 1e-6, 1e6)
        lambda2_new = np.clip(lambda2_new, 1e-6, 1e6)
    
        # Updated mixing proportions
        # pi_new[k] = average responsibility for component k
        pi_new = np.array([
            gamma1.sum() / n,
            gamma2.sum() / n
        ], dtype=float)
    
        # Store histories for diagnostics and convergence assessment
        loglik_history.append(loglik)
        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())
    
        # Parameter changes for convergence checking
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
    
        # Commit updates for next iteration
        lambda1 = lambda1_new
        lambda2 = lambda2_new
        pi = pi_new
    
        # %% Convergence check: stop if both log‑likelihood and parameters stabilize
        if it > 0:
            if (
                abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik
                and delta_lambda < tol_params
                and delta_pi < tol_params
            ):
                break
    
    # %% Final output: histories and last‑iteration parameter values
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
