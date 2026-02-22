# SAINT: Significance Analysis of INTeractome, Revisited

This repository contains a fully reproducible Python implementation of both **classical SAINT** and a **hierarchical EM-based SAINT model** for analyzing IP–MS interactome data. The package provides:

- A clean, modular pipeline for classical SAINT  
- A hierarchical mixture model solved via EM (not MCMC)  
- Automatic wide→long data handling  
- Deterministic gamma3-based scoring  
- Publication-quality diagnostics (volcano plots, density plots, staggered heatmaps)  
- A complete JupyterLab-friendly notebook workflow using # %% cell markers  

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

SAINT (Significance Analysis of Interactome Networks) is a probabilistic framework for identifying true protein–protein interactions from IP–MS data. This repository provides:

- **Classical SAINT**: A two-component mixture model with EM updates.  
- **Hierarchical SAINT**: A structured, multi-level mixture model with hierarchical priors on mixture weights and component means, solved via EM.

The hierarchical model is **deterministic**, **not MCMC**, and provides more stable gamma3 estimates through structured shrinkage.

---

## 2. Installation

Clone the repository:

```
git clone https://github.com/<eriklarsen4>/<Endo>.git
cd <Endo>
```

Install dependencies:

```
pip install -r requirements.txt
```

The package is pure Python and requires no compiled extensions.

---

## 3. Data Requirements

The pipeline accepts **wide-format** IP–MS data with:

- A `Protein` column containing gene names  
- A `MW` (molecular weight) column  
- Replicate columns for each bait, e.g.:

```
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
The pipeline automatically reshapes wide-format data into long-format per-bait matrices.

### Biological mapping rule  
Bait→protein mappings are **explicitly defined by the user** and treated as **gospel** by the pipeline.  
No inference or guessing is performed.

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

Both pipelines return a tidy results DataFrame with:

- gamma3 values  
- lambda parameters  
- mixture weights  
- replicate summaries  
- FDR estimates  

---

## 5. Hierarchical Model

### Developer-Level Explanation

The hierarchical SAINT model is a **three-component mixture model** with:

- Latent class assignments \( z_{i,b} \)  
- Component means \( \lambda_1, \lambda_2, \lambda_3 \)  
- Hierarchical priors on mixture weights \( \pi \)  
- Hierarchical priors on component means  
- Deterministic EM updates  

The algorithm alternates:

1. **E-step**: Compute responsibilities \( \gamma_{i,b,k} \)  
2. **M-step**: Update \( \lambda_k \), \( \pi_k \), and hyperparameters  
3. Repeat until convergence  

This is **not** MCMC.  
There are **no posterior draws**, **no chains**, and **no sampling**.

### Layperson Explanation

The hierarchical model groups proteins into three categories:

1. Background  
2. Weak binders  
3. True interactors  

It learns these categories by comparing replicate abundances across baits, borrowing strength across proteins to stabilize estimates.  
The gamma3 score is the probability that a protein belongs to the “true interactor” group.

---

## 6. Repository Structure

```
Endo/
│
├── saint/
│   ├── pipeline/
│   │   ├── classical_saint.py
│   │   ├── hierarchical_saint.py
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

## 7. Interpretation Guide

### gamma3  
Responsibility for the high-abundance component.

### FDR  
Computed deterministically from sorted gamma3 values.

### Volcano Plot  
- x-axis: inverted signal-to-noise ratio  
- y-axis: gamma3  
- Highlighted proteins: top interactors + TMEM184B  

### Staggered Heatmap  
Shows raw replicate counts for:

- TMEM184B  
- Top 25 interactors  
- Bottom 10 interactors  

---

## 8. Glossary

**gamma3** — Responsibility for the high-intensity component.  
**lambda1/2/3** — Component means.  
**pi** — Mixture weights.  
**E-step** — Computes expected latent assignments.  
**M-step** — Updates parameters.  
**Hierarchical prior** — Shrinkage structure across proteins.  
**Wide→long** — Automatic data reshaping.  

---

## 9. Appendix: Model Structure Summary

The hierarchical SAINT model is:

- A structured mixture model  
- Solved via EM  
- Deterministic  
- Stable across runs  
- More expressive than classical SAINT  
- Not a Bayesian sampler  

It provides a principled gamma3 score that integrates replicate structure, bait structure, and hierarchical shrinkage.
