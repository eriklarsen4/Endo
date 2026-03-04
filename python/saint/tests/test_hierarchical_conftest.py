# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:36:16 2026

@author: Erik
"""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def minimal_input_data():
    # Two baits, 5 proteins each, 2 replicates
    data = {
        "bait": ["B1"] * 10 + ["B2"] * 10,
        "prey": [f"P{i}" for i in range(1, 11)] * 2,
        "rep1": np.random.poisson(5, 20),
        "rep2": np.random.poisson(5, 20),
    }
    return pd.DataFrame(data)