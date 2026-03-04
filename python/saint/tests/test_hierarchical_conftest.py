# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:34:21 2026

@author: Erik
"""
# %% Import
import pytest
import pandas as pd

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline

# %% fixtures
@pytest.fixture
def minimal_input_data():
    return pd.DataFrame(
        {
            "Protein": ["P1", "P2"],
            "condA_B1_1": [5, 1],
            "condA_B1_2": [3, 0],
            "condA_CTRL_1": [2, 0],
            "condA_CTRL_2": [1, 1],
        }
    )


@pytest.fixture
def hierarchical_results(minimal_input_data):
    return run_hierarchical_pipeline(
        input_data=minimal_input_data,
        max_iter=5,
        seed=1,
        make_plots=False,
    )


@pytest.fixture
def hierarchical_metadata(hierarchical_results):
    return hierarchical_results["metadata"]


@pytest.fixture
def hierarchical_raw(hierarchical_results):
    return hierarchical_results["raw_outputs"]


@pytest.fixture
def hierarchical_df(hierarchical_results):
    return hierarchical_results["results_df"]