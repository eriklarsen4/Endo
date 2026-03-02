# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:32:12 2026

@author: Erik
"""

# %% Imports
import numpy as np
import pandas as pd
from saint.io.data_input import load_bait_data
from saint.model.hierarchical_tau_grid import run_tau_grid


# %% Hierarchical pipeline
def run_hierarchical_pipeline(
    input_data,
    bait_names=None,
    metadata=None,
    hyperparams=None,
    make_plots=False,
    plot_dir=None,
    max_iter=100,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=None,
    verbose=False,
):
    """
    Run the hierarchical SAINT pipeline for all baits using a per-bait tau grid
    search. This function loads the data, applies the hierarchical EM model with
    a two-stage tau grid for each bait, collects the best fitting results, and
    assembles a unified prey-level results table.

    Hyperparameters may be omitted entirely. If the user supplies no dictionary,
    or supplies a dictionary missing some entries, all missing hyperparameters
    are filled with internal defaults before any access occurs. Only the `tau`
    entry is modified during the grid search; all other hyperparameters remain
    unchanged unless explicitly provided by the user.

    Parameters
    
    input_data : dict or DataFrame
        Raw APMS data in the format expected by the data loader. The loader
        converts the wide-format APMS table into per-bait replicate matrices
        and constructs metadata including the per-bait protein ordering.

    hyperparams : dict, optional
        Optional hyperparameters for the hierarchical EM model. Users may supply
        overrides, but any missing entries are filled with internal defaults
        before the model is run. The `tau` entry is overridden during the grid
        search; all other entries remain unchanged unless explicitly provided.

    Returns
    
    dict
        Unified results object containing:
            - raw_outputs:
                em_results : dict mapping each bait to the best-fitting EM result
                tau_info   : dict mapping each bait to tau grid diagnostics
            - metadata : experiment metadata produced by the data loader
            - results_df : combined prey-level results table with one row per
              prey–bait pair, containing:
                  Protein, bait,
                  lambda1, lambda2, lambda3,
                  tau,
                  pi1, pi2, pi3,
                  gamma1, gamma2, gamma3

    Model parameters
    
    The model parameters lambda1, lambda2, lambda3, pi, and gamma are estimated
    by the EM algorithm and are not supplied by the user.

    Hyperparameters
    
    The hyperparameters alpha, a_k, b_k, and tau may be supplied in the
    hyperparameter dictionary but are not required. Any missing entries are
    filled with internal defaults before any access occurs. The `tau` entry is
    overridden by the grid search. All other hyperparameters remain unchanged
    unless explicitly provided.

    Tau grid tuning parameters
    
    The "coarse" and "refinement" grid ranges and resolutions define the search
    over tau. These are tuning parameters for the grid search procedure and are
    not model hyperparameters.
    """

    # %% Load data
    metadata_arg = metadata
    bait_list, X_by_bait, metadata = load_bait_data(input_data)
    loader_metadata = metadata

    # %% Storage
    all_results = {}
    all_tau_info = {}
    all_convergence_info = {}
    all_iteration_counts = {}

    # %% Default hyperparameters (must come BEFORE any access)
    if hyperparams is None:
        hyperparams = {}

    # Dirichlet prior (alpha = 2 for each component)
    hyperparams.setdefault("alpha1", 2.0)
    hyperparams.setdefault("alpha2", 2.0)
    hyperparams.setdefault("alpha3", 2.0)

    # Mixing proportions initialization (uniform)
    hyperparams.setdefault("pi_init", np.ones(3, dtype=float) / 3.0)

    # Tau grid (your actual coarse grid: logspace(-1, 1, 7))
    hyperparams.setdefault("tau_grid", np.logspace(-1.0, 1.0, num=7))

    # Tau default = geometric center of the grid = 10^0 = 1.0
    hyperparams.setdefault("tau", 1.0)

    # Now safe to access
    tau_grid = hyperparams["tau_grid"]

    # %% Per bait tau grid
    for bait in bait_list:
        X = X_by_bait[bait]

        tau_output = run_tau_grid(X, hyperparams, bait)

        all_results[bait] = tau_output["best_result"]
        all_tau_info[bait] = {
            "best_tau": tau_output["best_tau"],
            "tau_grid_results": tau_output["tau_grid_results"]
        }

        all_convergence_info[bait] = tau_output["convergence_info"]
        all_iteration_counts[bait] = tau_output["iteration_counts"]

    # %% Build results_df inline (no helper)
    rows = []

    for bait in bait_list:
        result = all_results[bait]
        proteins = metadata["proteins_by_bait"][bait]

        lambda1 = result["lambda1"]
        lambda2 = result["lambda2"]
        lambda3 = result["lambda3"]
        tau = result["tau"]
        pi = result["pi"]
        gamma = result["gamma"]

        df_bait = pd.DataFrame({
            "Protein": proteins,
            "bait": bait,
            "lambda1": lambda1,
            "lambda2": lambda2,
            "lambda3": lambda3,
            "tau": tau,
            "pi1": pi[0],
            "pi2": pi[1],
            "pi3": pi[2],
            "gamma1": gamma[:, 0],
            "gamma2": gamma[:, 1],
            "gamma3": gamma[:, 2],
        })

        rows.append(df_bait)

    results_df = pd.concat(rows, ignore_index=True)
    results_df = results_df.sort_values(["Protein", "gamma3"], ascending=[True, False])

    # %% Build the typed metadata

    from saint.pipeline.metadata_types import (
        HierarchicalMetadata,
        UserProvidedFields,
        InferredFields,
        PipelineDerivedFields,
    )

    # %% Construct typed metadata
    metadata_obj = HierarchicalMetadata(
        user_provided_fields=UserProvidedFields(
            biological_bait_names=metadata.get("biological_bait_names", {}),
            AN=metadata.get("AN", {}),
            MW=metadata.get("MW", {}),
            extra_fields={
                k: v for k, v in metadata.items()
                if k not in {"biological_bait_names", "AN", "MW"}
            },
        ),
        inferred_fields=InferredFields(
            baits=loader_metadata.get("baits", []),
            proteins_by_bait=loader_metadata.get("proteins_by_bait", {}),
            replicate_map=loader_metadata.get("replicate_map", {}),
            conditions=loader_metadata.get("conditions", {}),
            negative_controls_inferred=loader_metadata.get("negative_controls_inferred", []),
            extra_fields={
                k: v for k, v in loader_metadata.items()
                if k not in {
                    "baits",
                    "proteins_by_bait",
                    "replicate_map",
                    "conditions",
                    "negative_controls_inferred",
                }
            },
        ),
        pipeline_derived_fields=PipelineDerivedFields(
            bait_names=bait_names,
            hyperparameters=hyperparams,
            tau_grid=tau_grid,
            convergence=all_convergence_info,
            iteration_counts=all_iteration_counts,
        ),
    )
    
    # %% Plotting (optional)
    if make_plots:
         import matplotlib.pyplot as plt
         from saint.diagnostics.diagnostics_hierarchical import (
             make_hierarchical_plots,
             plot_gamma3_density,
             plot_multi_bait_summary,
         )
         from saint.diagnostics.diagnostics_tau_grid import diagnostics_tau_grid
    
         # Per-bait EM diagnostics
         for bait in bait_list:
             make_hierarchical_plots(
                 results_em=all_results[bait],
                 bait_name=bait,
             )
    
         # Tau-grid diagnostics (per bait)
         for bait in bait_list:
             diagnostics_tau_grid(
                 tau_info=all_tau_info[bait],
                 bait_name=bait,
             )
    
         # Gamma3 density across all baits
         plot_gamma3_density(results_df)
    
         # Multi-bait summary
         plot_multi_bait_summary(
             results_df=results_df,
             bait_list=bait_list,
         )
    
         plt.show()

    # %% Final unified output
    return {
        "raw_outputs": {
            "em_results": all_results,
            "tau_info": all_tau_info,
        },
        "metadata": metadata_obj,
        "results_df": results_df,
    }