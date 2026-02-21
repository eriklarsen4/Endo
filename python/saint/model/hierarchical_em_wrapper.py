# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 20:20:26 2026

@author: Erik
"""

# %% Imports

import numpy as np

from saint.model.hierarchical_responsibilities import compute_responsibilities
from saint.model.updates.lambda_updates import update_lambda_hierarchical
from saint.model.updates.tau_updates import update_tau
from saint.model.updates.pi_updates import update_pi
from saint.model.likelihood import compute_loglik


# %% Hierarchical EM wrapper

def run_em_hierarchical(
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
    Run hierarchical EM for one bait.

    lambda1, lambda2, lambda3 are component-specific Poisson rates
    for each prey (protein). These represent "background", "intermediate",
    and "signal" components, respectively, in a three-component mixture.

    pi is the vector of mixture proportions for the three components.
    gamma represents the posterior probability that each prey belongs
    to each component.

    alpha is the shape parameter of the Gamma prior on each lambda.
    tau is the shared rate parameter that controls global shrinkage
    of all lambda values.

    The E step computes gamma using the current lambda and pi values.
    The M step updates lambda1, lambda2, lambda3, tau, and pi using gamma.

    X is a numeric matrix of prey (protein) by bait.
    hyperparams is a dictionary containing initial lambda tau alpha and pi values.
    biological_bait is the protein symbol corresponding to the target bait.

    This wrapper is called once per bait by the multi bait pipeline.
    """

    if seed is not None:
        np.random.seed(seed)

    lambda1 = hyperparams["lambda1_init"].copy()
    lambda2 = hyperparams["lambda2_init"].copy()
    lambda3 = hyperparams["lambda3_init"].copy()
    tau = hyperparams["tau_init"].copy()
    pi = hyperparams["pi_init"].copy()
    alpha = hyperparams["alpha"].copy()

    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    lambda3_history = []
    tau_history = []
    pi_history = []
    gamma_history = []

    prev_loglik = None

    for iteration in range(max_iter):

        gamma = compute_responsibilities(
            X,
            lambda1,
            lambda2,
            lambda3,
            pi
        )

        lambda1_new, lambda2_new, lambda3_new = update_lambda_hierarchical(
            X,
            gamma,
            alpha,
            tau
        )

        tau_new = update_tau(
            X,
            gamma,
            lambda1_new,
            lambda2_new,
            lambda3_new
        )

        pi_new = update_pi(gamma)

        loglik = compute_loglik(
            X,
            lambda1_new,
            lambda2_new,
            lambda3_new,
            pi_new
        )
        loglik_history.append(loglik)

        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        lambda3_history.append(lambda3_new.copy())
        tau_history.append(tau_new.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())

        delta_lambda1 = np.max(np.abs(lambda1_new - lambda1))
        delta_lambda2 = np.max(np.abs(lambda2_new - lambda2))
        delta_lambda3 = np.max(np.abs(lambda3_new - lambda3))
        delta_tau = np.max(np.abs(tau_new - tau))
        delta_pi = np.max(np.abs(pi_new - pi))

        param_change = max(
            delta_lambda1,
            delta_lambda2,
            delta_lambda3,
            delta_tau,
            delta_pi
        )

        if verbose:
            print(
                f"Iter {iteration}  loglik {loglik:.6f}  param_change {param_change:.6e}"
            )

        if prev_loglik is not None:
            if (
                abs(loglik - prev_loglik) < tol_loglik
                and param_change < tol_params
            ):
                break

        lambda1 = lambda1_new
        lambda2 = lambda2_new
        lambda3 = lambda3_new
        tau = tau_new
        pi = pi_new
        prev_loglik = loglik

    results = {
        "lambda1": lambda1,
        "lambda2": lambda2,
        "lambda3": lambda3,
        "tau": tau,
        "pi": pi,
        "gamma": gamma,
        "loglik_history": np.array(loglik_history),
        "lambda1_history": lambda1_history,
        "lambda2_history": lambda2_history,
        "lambda3_history": lambda3_history,
        "tau_history": tau_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history,
        "biological_bait": biological_bait
    }

    return results



