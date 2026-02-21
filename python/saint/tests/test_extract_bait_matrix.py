# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:20:16 2026

@author: Erik
"""

# %% Imports

import pandas as pd

from saint.io.data_input import extract_bait_matrix


# %% Fixtures

def make_input():
    """
    Create a small input table with two baits.
    """

    data = {
        "Protein": ["P1", "P2", "P3", "P1", "P2", "P3"],
        "Bait": ["B1", "B1", "B1", "B2", "B2", "B2"],
        "Count": [10, 2, 0, 5, 1, 0]
    }

    return pd.DataFrame(data)


# %% Tests

def test_extract_bait_matrix_shapes():
    """
    Test that extract_bait_matrix returns the expected shapes.
    """

    df = make_input()
    X, metadata, column_labels, mapping = extract_bait_matrix(df, "B1")

    assert X.shape == (3, 1)
    assert metadata.shape[0] == 3
    assert len(column_labels) == 1


def test_extract_bait_matrix_filters_correct_bait():
    """
    Test that extract_bait_matrix selects only the requested bait.
    """

    df = make_input()
    X, metadata, column_labels, mapping = extract_bait_matrix(df, "B2")

    assert metadata["Protein"].tolist() == ["P1", "P2", "P3"]
    assert X[:, 0].tolist() == [5, 1, 0]