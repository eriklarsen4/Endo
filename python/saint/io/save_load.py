# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 20:13:57 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd
from pathlib import Path

# %% Save

def save_results(results, output_path):
    """
    Save pipeline results to disk.

    Parameters
    results : dict
        Contains output_df and raw_outputs
    output_path : string or Path
    """

    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    output_df = results["output_df"]
    raw_outputs = results["raw_outputs"]

    output_df.to_csv(output_path / "saint_results.csv", index=False)

    np.save(output_path / "raw_outputs.npy", raw_outputs, allow_pickle=True)

# %% Load

def load_results(output_path):
    """
    Load pipeline results from disk.

    Parameters
    output_path : string or Path

    Returns
    results : dict
        Contains output_df and raw_outputs
    """

    output_path = Path(output_path)

    output_df = pd.read_csv(output_path / "saint_results.csv")

    raw_outputs = np.load(output_path / "raw_outputs.npy", allow_pickle=True).item()

    return {
        "output_df": output_df,
        "raw_outputs": raw_outputs
    }