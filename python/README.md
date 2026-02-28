# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 08:07:02 2026

@author: Erik
"""

# SAINT: Significance Analysis of INTeractome, Revisited

This repository contains a fully reproducible Python implementation of both **classical SAINT** and a **hierarchical EM-based SAINT model** for analyzing IP–MS interactome data. The package provides:

- A clean, modular pipeline for classical SAINT  
- A hierarchical three-component mixture model solved via deterministic EM  
- Automatic wide→long data handling  
- Structured shrinkage via empirical Bayes  
- Deterministic gamma3-based scoring  
- Publication-quality diagnostics (volcano plots, density plots, staggered heatmaps)  
- A JupyterLab-friendly workflow using `# %%` cell markers  

---

## Table of Contents

1. Overview  
2. Installation  
3. Data Requirements  
4. Running the Pipelines  
5. Hierarchical Model (Developer + Layperson Explanations)  
6. Notebook Workflow  
7. Repository Structure  
8. Interpretation Guide  
9. Glossary  
10. Appendix: Model Structure Summary  

---

## 1. Overview

SAINT (Significance Analysis of INTeractome) is a probabilistic framework for identifying true protein–protein interactions from IP–MS data. This repository implements:

- **Classical SAINT** — a two-component Poisson mixture model with EM updates  
- **Hierarchical SAINT** — a structured three-component mixture model with hierarchical priors on mixture weights and component means, solved via EM  

The hierarchical model is **deterministic**.  
There is **no MCMC**, **no sampling**, and **no posterior chains**.  
All updates are closed-form and stable across runs.

---

## 2. Installation

Clone the repository:

```bash
git clone https://github.com/eriklarsen4/Endo.git
cd Endo
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Pure Python. No compiled extensions.

---

## 3. Data Requirements

The pipeline accepts **wide-format** IP–MS data with:

- `Protein` column  
- `MW` (molecular weight) column  
- replicate columns for each bait, e.g.:

```text
Control_BirAV5_1
Control_BirAV5_2
Control_GFPmyc_1
Control_GFPmyc_2
Treat_TMEMV5_1
Treat_TMEMV5_2
Treat_TMEMmyc_1
Treat_TMEMmyc_2
```

### Automatic wide→long conversion  
The pipeline automatically reshapes wide-format data into per-bait replicate matrices.

### Explicit bait mapping  
Users provide the bait→protein mapping.  
The pipeline **never infers** or guesses bait identity.

---

## 4. Running the Pipelines

### Classical SAINT

```python
from saint.pipeline.classical_saint import run_classical_pipeline

results_classical = run_classical_pipeline(
    input_data=data,
    bait_names=["TMEMV5", "TMEMmyc"],
    metadata=metadata,
    make_plots=True
)
```

### Hierarchical SAINT (EM-based)

```python
from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline

results = run_hierarchical_pipeline(
    input_data=data,
    bait_names=["TMEMV5", "TMEMmyc"],
    metadata=metadata,
    make_plots=True
)
```

Both pipelines return a tidy DataFrame containing:

- gamma3 values  
- λ₁/λ₂/λ₃ estimates  
- mixture weights  
- replicate summaries  
- FDR estimates  

---

## 5. Hierarchical Model

### Developer-Level Explanation

The hierarchical SAINT model is a **three-component Poisson mixture** with:

- latent class assignments  
- component means λ₁, λ₂, λ₃  
- Dirichlet prior on mixture weights π  
- Gamma(aₖ, bₖ) priors on component means  
- empirical Bayes updates for all hyperparameters  
- deterministic EM updates  

Component semantics:

- **λ₁** — background  
- **λ₂** — contaminant / weak binder  
- **λ₃** — true signal  

The EM algorithm alternates:

1. **E-step**  
   Compute responsibilities γ using stabilized log-sum-exp.

2. **Empirical Bayes updates**  
   - Dirichlet parameters updated by moment matching on γ  
   - Gamma(aₖ, bₖ) hyperparameters updated by moment matching on λₖ  
     - mean = a / b  
     - variance = a / b²  

3. **M-step**  
   Closed-form updates for λ₁, λ₂, λ₃ and π (method of moments).

This is a **deterministic hierarchical EM**, not a Bayesian sampler.

### Layperson Explanation

The hierarchical model groups proteins into:

1. Background  
2. Weak binders  
3. True interactors  

It learns these groups by comparing replicate counts across proteins and baits, borrowing strength across proteins to stabilize estimates.  
The gamma3 score is the probability that a protein belongs to the “true interactor” group.

---

## 6. Notebook Workflow

All pipeline modules use `# %%` cell markers for:

- clean JupyterLab execution  
- stepwise debugging  
- reproducible workflows  

---

## 7. Repository Structure

```text
Endo/
│
├── saint/
│   ├── pipeline/
│   │   ├── classical_saint.py
│   │   ├── hierarchical_saint.py
│   │   └── ...
│   ├── model/
│   │   ├── classical_em_wrapper.py
│   │   ├── hierarchical_em_wrapper.py
│   │   └── ...
│   ├── diagnostics/
│   ├── utils/
│   └── __init__.py
│
├── data/
│   └── IPMS_counts.rda
│
├── notebooks/
│   └── proteomics.ipynb
│
├── README.md
└── requirements.txt
```

---

## 8. Interpretation Guide

### gamma3  
Responsibility for the **signal** component (λ₃).

### FDR  
Computed deterministically from sorted gamma3 values.

### Volcano Plot  
- x-axis: inverted signal-to-noise  
- y-axis: gamma3  
- highlights: top interactors + TMEM184B  

### Staggered Heatmap  
Shows replicate-level counts for:

- TMEM184B  
- top 25 interactors  
- bottom 10 interactors  

---

## 9. Glossary

- **gamma3** — probability of belonging to the signal component  
- **λ₁/λ₂/λ₃** — background, contaminant, signal rates  
- **π** — mixture weights  
- **E-step** — expected latent assignments  
- **M-step** — parameter updates  
- **hierarchical prior** — shrinkage across proteins  
- **moment matching** — closed-form update for Gamma priors  
- **wide→long** — automatic data reshaping  

---

## 10. Appendix: Model Structure Summary

The hierarchical SAINT model is:

- a structured three-component mixture  
- solved via deterministic EM  
- equipped with empirical Bayes shrinkage  
- stable across runs  
- more expressive than classical SAINT  
- not a Bayesian sampler  

It produces a principled gamma3 score that integrates replicate structure, bait structure, and hierarchical shrinkage.