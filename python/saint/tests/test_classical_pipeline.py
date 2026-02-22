# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:07:19 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd

from saint.pipeline.classical_saint import run_classical_pipeline


# %% Dummy EM wrapper

class DummyEM:
    def __init__(self, n):
        self.n = n

    def run(self, X, init_params, biological_bait, max_iter, tol_loglik, tol_params, seed, verbose):
        return {
            "loglik_history": [0.0, 1.0],
            "lambda1_history": [],
            "lambda2_history": [],
            "pi_history": [],
            "gamma_history": [],
            "lambda1": np.ones(self.n),
            "lambda2": np.ones(self.n) * 2,
            "pi": np.tile(np.array([0.3, 0.7]), (self.n, 1)),
            "gamma": np.tile(np.array([0.1, 0.9]), (self.n, 1))
        }


# %% Test classical pipeline

def test_classical_pipeline(monkeypatch):
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condX_BaitA_1": [5, 1],
        "condX_BaitA_2": [7, 0],
        "condX_Ctrl1_1": [2, 0],
        "condX_Ctrl1_2": [1, 1]
    })

    metadata = {
        "biological_bait_names": {
            "BaitA": "BioA"
        },
        "AN": {"BioA": "A1"},
        "MW": {"BioA": 10000}
    }

    dummy = DummyEM(n=2)
    monkeypatch.setattr(
        "saint.model.classical_em_wrapper.run_em",
        dummy.run
    )

    results = run_classical_pipeline(
        df,
        bait_names=["BaitA"],
        metadata=metadata
    )

    results_df = results["results_df"]

    assert "lambda1" in results_df.columns
    assert "gamma2" in results_df.columns

    assert results_df.iloc[0]["Protein"] == "P1"
    assert results_df.iloc[0]["gamma2"] == 0.9

    assert "controls_used" in results["metadata"]
    assert "Ctrl1" in results["metadata"]["controls_used"]


