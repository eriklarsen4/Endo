# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:12:39 2026

@author: Erik
"""

# Future Tests Roadmap

These tests are intentionally deferred until the hierarchical and classical
pipelines are fully stabilized and integrated with real EM behavior.

## 1. Missing replicate handling
Validate that extract_bait_matrix correctly handles:
- missing rep2 when rep1 exists
- uneven replicate counts across baits
- non-sequential replicate labels

## 2. Multi-bait sorting stability
Validate that sorting by gamma3 (hierarchical) and gamma2 (classical) is stable
within each Protein group when multiple baits exist.

## 3. Metadata integrity
Validate that metadata_out contains:
- bait_names
- biological_bait_names
- AN
- MW
- controls_used
and that no unexpected keys appear.

## 4. Endogenous datasets
Validate that extract_bait_matrix and both pipelines behave correctly when
bait names equal biological bait names and no tags are present.

## 5. EM wrapper error handling
Validate that both pipelines raise clean exceptions when:
- EM returns NaN parameters
- EM returns mismatched array lengths
- EM fails to converge

## 6. Controls edge cases
Validate behavior when:
- no controls exist
- all baits are controls
- controls have missing replicates
