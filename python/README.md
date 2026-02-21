# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 16:34:29 2026

@author: Erik
"""

<p align="center">

<a href="https://github.com/eriklarsen4/Endo/actions">
<img src="https://img.shields.io/github/actions/workflow/status/eriklarsen4/Endo/tests.yml?label=tests" alt="Tests">
</a>

<a href="https://codecov.io/gh/eriklarsen4/Endo">
<img src="https://img.shields.io/codecov/c/github/eriklarsen4/Endo" alt="Coverage">
</a>

<img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python Version">

<img src="https://img.shields.io/badge/license-MIT-green" alt="License">

</p>

# SAINT

SAINT is a modular, explicit, and reproducible implementation of classical and hierarchical SAINT models.  
The package follows a strict layered architecture:

```
pipeline → wrapper → EM engine → updates and responsibilities → diagnostics
```

All components are literal, minimal, and future proof.  
There is no hidden state and no implicit inference beyond what is documented.

---

## Table of contents

- [Overview](#overview)
- [Model summary](#model-summary)
- [Workflow diagram](#workflow-diagram)
- [Input data](#input-data)
- [Modes of operation](#modes-of-operation)
  - [Default mode](#default-mode)
  - [Override mode](#override-mode)
- [Classical pipeline](#classical-pipeline)
- [Hierarchical pipeline](#hierarchical-pipeline)
- [Diagnostics](#diagnostics)
- [Testing](#testing)
- [Reproducibility](#reproducibility)
- [Usage example](#usage-example)
- [How to interpret SAINT outputs](#how-to-interpret-saint-outputs)
- [Glossary](#glossary)
- [Appendix A — Developer level EM explanation](#appendix-a--developer-level-em-explanation)
- [Appendix B — Layperson biological explanation](#appendix-b--layperson-biological-explanation)
- [Appendix C — Workflow diagram](#appendix-c--workflow-diagram)
- [License](#license)

---

## Overview

SAINT fits Poisson mixture models to quantitative interaction data for each bait.  
The package provides:

- explicit data extraction through `extract_bait_matrix`
- classical and hierarchical pipelines (`run_classical_pipeline`, `run_hierarchical_pipeline`)
- EM orchestration (`run_em`, `run_em_hierarchical`)
- pure mathematical updates in `saint/model/updates`
- diagnostics that return matplotlib Figure objects for interactive inspection

Two modes of operation are supported (see [Modes of operation](#modes-of-operation)):

- **Default mode**
- **Override mode**

Mode selection is determined by the presence of the metadata key `bait_to_biological_name`.

---

## Model summary

SAINT uses mixture models with prey specific Poisson rates:

- **lambda1**: background rate  
- **lambda2**: signal rate (classical)  
- **lambda3**: intermediate component (hierarchical only)  
- **pi**: mixture proportions  
- **gamma**: posterior component probabilities  
- **alpha**: Gamma prior shape parameter (hierarchical)  
- **tau**: Gamma prior rate parameter controlling global shrinkage (hierarchical)

The EM algorithm alternates between:

- **E step**: compute `gamma` (responsibilities)  
- **M step**: update `lambda`, `pi`, and `tau` (hierarchical only)

For more detail, see [Appendix A — Developer level EM explanation](#appendix-a--developer-level-em-explanation).  
For a biological interpretation, see [Appendix B — Layperson biological explanation](#appendix-b--layperson-biological-explanation).

---

## Workflow diagram

```
Input DataFrame
    ↓
extract_bait_matrix
    ↓
pipeline (classical or hierarchical)
    ↓
wrapper (run_em or run_em_hierarchical)
    ↓
EM engine
    - responsibilities (gamma)
    - lambda updates
    - pi updates
    - tau updates (hierarchical)
    - log likelihood
    ↓
diagnostics (matplotlib Figures)
    ↓
output_df + histories + metadata
```

A more detailed narrative is given in  
[Appendix C — Workflow diagram](#appendix-c--workflow-diagram).

---

## Input data

The primary input is a pandas DataFrame or a CSV path.  
Required columns:

- `Protein`
- `Bait`
- one or more `Count_*` columns

If the input is in wide format (Protein × bait columns), the pipeline automatically converts it to long format using `extract_bait_matrix`.

---

## Modes of operation

SAINT supports two modes of operation that control how bait names and biological protein symbols are handled.

### Default mode

Default mode is active when:

- `metadata` is `None`, or  
- `metadata` does not contain `bait_to_biological_name`

In this mode:

- bait names are treated as biological protein symbols  
- no user provided mapping is required  
- behavior matches original SAINT conventions  

### Override mode

Override mode is active when:

- `metadata` contains `bait_to_biological_name`

In this mode:

- the user provides a mapping from bait names to biological protein symbols  
- this mapping is validated for all baits  
- positive control prey is the biological protein symbol  
- negative control inference still uses column names  
- override mode is required when bait names differ from protein symbols  

---

## Classical pipeline

Run through `run_classical_pipeline`.

The classical pipeline performs:

- mode detection (default or override)
- metadata validation
- bait specific extraction
- hyperparameter initialization
- EM fitting through `run_em`
- diagnostics via `make_classical_plots`
- output assembly

The output is a dictionary keyed by bait name containing:

- `output_df`: final lambda1, lambda2, pi, and gamma values  
- `merged_df`: input data merged with `output_df`  
- `raw_outputs`: histories of lambda, pi, gamma, and log likelihood  
- `metadata`: bait and biological bait names  

See [How to interpret SAINT outputs](#how-to-interpret-saint-outputs) for guidance on reading `output_df`.

---

## Hierarchical pipeline

Run through `run_hierarchical_pipeline`.

This pipeline mirrors the classical pipeline but includes:

- a third lambda component (`lambda3`)
- a prey specific shrinkage parameter (`tau`)
- Gamma priors defined by `alpha` and `tau`
- hierarchical EM fitting through `run_em_hierarchical`

The output structure matches the classical pipeline, with additional columns for `lambda3` and `tau`.

---

## Diagnostics

Diagnostics live in:

- `diagnostics_classical.py`
- `diagnostics_hierarchical.py`

Each diagnostic function returns a dictionary of matplotlib Figure objects.  
Figures are not shown automatically, supporting interactive inspection in Spyder.

Classical diagnostics include:

- log likelihood history  
- lambda1 history  
- lambda2 history  
- pi history  
- gamma history  

Hierarchical diagnostics additionally include:

- lambda3 history  
- tau history  

These plots are useful when assessing convergence and stability of the EM algorithm.

---

## Testing

Tests live in `saint/tests`.  
Each test file has a single responsibility and mirrors the structure of its module.

- Coverage is configured through `.coveragerc`  
- Pytest behavior is configured through `pytest.ini`  

All directories contain an `__init__.py` file to ensure stable imports.

---

## Reproducibility

SAINT is designed for reproducible scientific workflows:

- all randomness is controlled through the `seed` parameter  
- all functions are deterministic given their inputs  
- no global state is used  
- pipelines are explicit and side effect free  

---

## Usage example

```python
from saint.pipeline.classical_saint import run_classical_pipeline

results = run_classical_pipeline(
    input_data="input.csv",
    bait_names=["BaitA", "BaitB"],
    metadata=None,
    max_iter=100
)
```

The output is a dictionary keyed by bait name.  
Each entry contains `output_df`, `merged_df`, `raw_outputs`, and `metadata`.

---

## How to interpret SAINT outputs

This section explains how to read the main outputs from the pipelines.  
See also [Glossary](#glossary) and [Appendix B](#appendix-b--layperson-biological-explanation).

For each bait, `output_df` contains one row per prey (protein) and columns such as:

- `lambda1`: estimated background rate for the prey  
- `lambda2`: estimated signal rate (classical)  
- `lambda3`: intermediate component rate (hierarchical only)  
- `pi1`, `pi2`, `pi3`: mixture proportions for each component  
- `gamma1`, `gamma2`, `gamma3`: posterior probabilities that the prey belongs to each component  

Typical interpretation:

- a prey with high `gamma2` (classical) or high `gamma3` (hierarchical) is more likely to be a true interactor  
- `lambda` values describe expected abundance under each component  
- `pi` values describe how common each component is across all preys  

The histories in `raw_outputs` (lambda, pi, gamma, tau, log likelihood) can be used to:

- check convergence  
- inspect stability of parameter estimates  
- diagnose potential issues with initialization or data quality  

---

## Glossary

- **Bait**: the target protein used to pull down interacting proteins  
- **Prey**: a protein detected in the experiment, potentially interacting with the bait  
- **X**: count matrix of prey by experiment for a given bait  
- **lambda1, lambda2, lambda3**: prey specific Poisson rates for background, intermediate, and signal components  
- **pi**: mixture proportions for the components  
- **gamma**: posterior probabilities that each prey belongs to each component  
- **alpha**: shape parameter of the Gamma prior on lambda (hierarchical)  
- **tau**: rate parameter of the Gamma prior controlling shrinkage (hierarchical)  
- **E step**: expectation step computing gamma  
- **M step**: maximization step updating lambda, pi, and tau  

For mathematical details, see [Appendix A](#appendix-a--developer-level-em-explanation).  
For biological intuition, see [Appendix B](#appendix-b--layperson-biological-explanation).

---

## Appendix A — Developer level EM explanation

This appendix summarizes the mathematical structure of the classical and hierarchical SAINT models.

Model structure:  
For each prey i and bait j, the observed count X_ij is modeled as a Poisson random variable with a component specific rate. Classical SAINT uses two components (background and signal). Hierarchical SAINT uses three components (background, intermediate, signal).

Lambda parameters:  
lambda1_i, lambda2_i, lambda3_i are prey specific Poisson rates. These represent the expected count for prey i under each mixture component.

Mixture proportions:  
pi_k are the mixture weights for components k. They satisfy pi_k ≥ 0 and sum_k pi_k = 1.

Responsibilities:  
gamma_ik are the posterior probabilities that prey i belongs to component k. They are computed in the E step using Bayes rule.

E step:  
gamma_ik ∝ pi_k * Poisson(X_i· | lambda_k_i)

M step (classical):  
lambda1_i = sum_j X_ij / sum_j gamma_i1  
lambda2_i = sum_j X_ij / sum_j gamma_i2  
pi_k = (1 / n_preys) * sum_i gamma_ik

M step (hierarchical):  
lambda_k_i = (sum_j X_ij + (alpha_k − 1)) / (n_conditions * gamma_ik + tau_i)

Tau update:  
tau_i = weighted_expected_counts_i / sum_k gamma_ik

Log likelihood:  
The observed data log likelihood is computed using the mixture of Poisson components.  
Convergence is determined by changes in log likelihood and parameter norms.

---

## Appendix B — Layperson biological explanation

This appendix explains the SAINT model in biological terms.

Goal:  
SAINT identifies which proteins (preys) are likely to interact with a target protein (bait) based on quantitative mass spectrometry data.

Counts:  
For each prey, we measure how many times it appears in experiments involving a specific bait. Higher counts suggest stronger or more consistent interactions.

Mixture model:  
SAINT assumes that each prey belongs to one of several biological categories:  
background (noise or non specific binding), intermediate (hierarchical only), and signal (likely true interactors).

Each category has a characteristic expected count level, represented by lambda values.

Posterior probabilities:  
For each prey, SAINT computes the probability that it belongs to each category. These probabilities are called gamma values.

Mixture proportions:  
The pi values describe how common each category is across all preys.

Hierarchical shrinkage:  
In the hierarchical model, tau controls how strongly the lambda values are pulled toward a common baseline. This prevents unstable or extreme estimates when data are sparse.

Interpretation:  
A prey with high gamma for the signal component is likely to be a true interactor.  
The lambda values describe expected abundance, and the pi values describe the overall distribution of categories.

---

## Appendix C — Workflow diagram

This appendix shows the full SAINT workflow from input to output.

Input:  
A DataFrame containing Protein, Bait, and Count_* columns.

Step 1: Data extraction  
extract_bait_matrix converts wide format to long format if needed and extracts bait specific count matrices.

Step 2: Pipeline  
The classical or hierarchical pipeline:
- validates metadata  
- initializes hyperparameters  
- calls the EM wrapper  
- assembles outputs  
- generates diagnostics  

Step 3: EM wrapper  
The wrapper runs the EM algorithm:
- E step: compute responsibilities (gamma)  
- M step: update lambda, pi, and tau (hierarchical only)  
- compute log likelihood  
- check convergence  
- store histories  

Step 4: Diagnostics  
Diagnostic functions return matplotlib Figure objects showing:
- log likelihood progression  
- lambda histories  
- pi histories  
- gamma histories  
- tau and lambda3 histories (hierarchical)  

Step 5: Output assembly  
The pipeline returns:
- output_df with final lambda, pi, gamma, and tau (hierarchical)  
- merged_df combining input and output  
- raw histories  
- metadata  

---

## License

MIT License.
