# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:25:09 2026

@author: Erik
"""

# %% Import
import matplotlib.figure as mpl_fig
import matplotlib.axes as mpl_ax

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline
from saint.diagnostics.diagnostics_hierarchical import (
    make_hierarchical_plots,
    plot_gamma3_density,
)

# %% test the diagnostics signature
def test_hierarchical_diagnostics(minimal_input_data):
    out = run_hierarchical_pipeline(
        input_data=minimal_input_data,
        make_plots=False,
        verbose=False,
    )

    for bait, em in out["raw_outputs"]["em_results"].items():
        figs = make_hierarchical_plots(em, bait)
        assert isinstance(figs, dict)
        for f in figs.values():
            assert isinstance(f, mpl_fig.Figure)

# %% test the density signature
def test_global_gamma3_density(minimal_input_data):
    out = run_hierarchical_pipeline(
        input_data=minimal_input_data,
        make_plots=False,
        verbose=False,
    )

    df = out["results_df"]
    ax = plot_gamma3_density(df)
    assert isinstance(ax, mpl_ax.Axes)