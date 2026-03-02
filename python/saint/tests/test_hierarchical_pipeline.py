# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:08:51 2026

@author: Erik
"""

# %% Imports
import pandas as pd
from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline
from saint.pipeline.metadata_types import (
    HierarchicalMetadata,
    UserProvidedFields,
    InferredFields,
    PipelineDerivedFields,
)


# %% Tests

def test_hierarchical_pipeline_runs_and_returns_expected_structure():
    # Minimal synthetic wide-format input
    input_data = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condA_B1_1": [5, 1],
        "condA_B1_2": [3, 0],
        "condA_CTRL_1": [2, 0],
        "condA_CTRL_2": [1, 1],
    })

    bait_names = ["B1"]

    # User metadata (not required to appear in output)
    metadata = {
        "biological_bait_names": {"B1": "BAIT1"},
        "AN": {"B1": "AN1"},
        "MW": {"B1": 50.0},
    }

    results = run_hierarchical_pipeline(
        input_data,
        bait_names,
        metadata,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    # Top-level keys
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata is a typed object
    meta = results["metadata"]
    assert isinstance(meta, HierarchicalMetadata)
    assert isinstance(meta.user_provided_fields, UserProvidedFields)
    assert isinstance(meta.inferred_fields, InferredFields)
    assert isinstance(meta.pipeline_derived_fields, PipelineDerivedFields)

    # Inferred fields must contain baits and proteins_by_bait
    assert len(meta.inferred_fields.baits) > 0
    assert isinstance(meta.inferred_fields.proteins_by_bait, dict)

    # Pipeline-derived fields must contain bait_names
    assert meta.pipeline_derived_fields.bait_names == bait_names

    # Raw outputs structure
    raw = results["raw_outputs"]
    assert "em_results" in raw
    assert "tau_info" in raw

    # Results dataframe
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)

    expected_cols = {
        "Protein", "bait",
        "lambda1", "lambda2", "lambda3",
        "tau",
        "pi1", "pi2", "pi3",
        "gamma1", "gamma2", "gamma3",
    }
    assert expected_cols.issubset(df.columns)

    # Sorting check: Protein ascending, gamma3 descending
    df_sorted = df.sort_values(
        by=["Protein", "gamma3"],
        ascending=[True, False],
        ignore_index=True,
    )
    assert df.equals(df_sorted)