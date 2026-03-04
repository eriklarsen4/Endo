# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:25:08 2026

@author: Erik
"""

# %% Import
import os
import tempfile

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline

# %% Test end to end works
def test_pipeline_end_to_end_with_plots(minimal_input_data):
    with tempfile.TemporaryDirectory() as tmp:
        out = run_hierarchical_pipeline(
            input_data=minimal_input_data,
            make_plots=True,
            show_plots=False,
            save_plots=True,
            plot_dir=tmp,
            verbose=False,
        )

        # Check that at least one hierarchical plot was saved
        saved = os.listdir(tmp)
        assert any("hier" in f for f in saved)

        # Check tau-grid plot saved
        assert any("tau_grid" in f for f in saved)

        # Check global gamma3 density saved
        assert any("global_gamma3_density" in f for f in saved)