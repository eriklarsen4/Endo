# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:08:51 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline


# %% Dummy EM wrapper

class DummyEM:
    def __init__(self, n):
        self.n = n

    def run(self, X, hyperparams, biological_bait, max_iter, tol_loglik, tol_params, seed, verbose):
        return {
            "loglik_history": [0.0, 1.0],
            "lambda1_history": [],
            "lambda2_history": [],
            "lambda3_history": [],
            "tau_history": [],
            "pi_history": [],
            "gamma_history": [],
            "lambda1": np.ones(self.n),
            "lambda2": np.ones(self.n) * 2,
            "lambda3": np.ones(self.n) * 3,
            "tau": np.ones(self.n) * 0.5,
            "pi": np.tile(np.array([0.1, 0.2, 0.7]), (self.n, 1)),
            "gamma": np.tile(np.array([0.05, 0.10, 0.85]), (self.n, 1))
        }


# %% Test hierarchical pipeline

def test_hierarchical_pipeline(monkeypatch):
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condX_BaitA_1": [5, 1],
        "condX_BaitA_2": [7, 0],
        "condX_BaitB_1": [3, 0],
        "condX_BaitB_2": [4, 1],
        "condX_Ctrl1_1": [2, 0],
        "condX_Ctrl1_2": [1, 1]
    })

    metadata = {
        "biological_bait_names": {
            "BaitA": "BioA",
            "BaitB": "BioB"
        },
        "AN": {"BioA": "A1", "BioB": "B1"},
        "MW": {"BioA": 10000, "BioB": 20000}
    }

    dummy = DummyEM(n=2)
    monkeypatch.setattr(
        "saint.model.hierarchical_em_wrapper.run_em",
        dummy.run
    )

    results = run_hierarchical_pipeline(
        df,
        bait_names=["BaitA", "BaitB"],
        metadata=metadata
    )

    results_df = results["results_df"]

    assert "rep1" in results_df.columns
    assert "rep2" in results_df.columns
    assert "gamma3" in results_df.columns

    assert results_df.iloc[0]["Protein"] == "P1"
    assert results_df.iloc[0]["gamma3"] == 0.85

    assert "controls_used" in results["metadata"]
    assert "Ctrl1" in results["metadata"]["controls_used"]
