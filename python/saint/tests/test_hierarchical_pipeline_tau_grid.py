# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:25:09 2026

@author: Erik
"""

# %% Import
from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline
from saint.diagnostics.diagnostics_tau_grid import diagnostics_tau_grid

# %% test
def test_tau_grid_and_best_tau(minimal_input_data):
    out = run_hierarchical_pipeline(
        input_data=minimal_input_data,
        make_plots=False,
        verbose=False,
    )

    tau_info = out["raw_outputs"]["tau_info"]
    assert isinstance(tau_info, dict)

    for bait, tg in tau_info.items():
        assert "taus" in tg
        assert "logliks" in tg
        assert len(tg["taus"]) == len(tg["logliks"])

        # Check best tau selection
        best_tau = out["raw_outputs"]["em_results"][bait]["tau"]
        assert best_tau in tg["taus"]

        # Check diagnostics interface
        diag = diagnostics_tau_grid(tg, bait, best_tau)
        assert "summary_df" in diag
        assert "figure" in diag