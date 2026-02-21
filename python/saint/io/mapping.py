# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 19:56:43 2026

@author: Erik
"""

import re

def build_mapping(column_labels):
    """
    Parse spectral column labels of the form <condition>_<bait>_<rep>.
    Ignore all non-matching columns (metadata).
    """

    pattern = re.compile(r"([^_]+)_([^_]+)_([^_]+)")

    all_columns = []
    by_condition = {}
    by_bait = {}
    by_pair = {}

    for col in column_labels:
        m = pattern.fullmatch(col)
        if m is None:
            # metadata or non-spectral column → skip
            continue

        condition, bait, rep = m.groups()

        all_columns.append(col)
        by_condition.setdefault(condition, []).append(col)
        by_bait.setdefault(bait, []).append(col)
        by_pair.setdefault((condition, bait), []).append(col)

    return {
        "all_columns": all_columns,
        "by_condition": by_condition,
        "by_bait": by_bait,
        "by_pair": by_pair
    }