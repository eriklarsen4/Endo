# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 21:09:54 2026

@author: Erik
"""

# %%  Import

import rdata
import os, sys
import pandas as pd

# Force working directory to the directory containing this script
# Spyder sets sys.argv[0] to the script path
script_path = os.path.abspath(sys.argv[0])
project_root = os.path.dirname(script_path)

os.chdir(project_root)
sys.path.insert(0, project_root)


# CHANGE THIS to the folder that contains the "saint" directory
#project_root = Path(r"C:\Users\Erik\Desktop\Programming\Python\Interactome Project")
#sys.path.insert(0, str(project_root))

# get the proteomics data
data = rdata.read_rda('C:/Users/Erik/Desktop/Programming/R/Bio/Endo/Rpkg/data/IPMS_counts.rda')

# %% 
# Data cleaning


# extract column names
names = list(pd.DataFrame(data['IPMS_counts'][:1]).iloc[0])

# convert the data from dictionary to pandas dataframe
data = pd.DataFrame(data['IPMS_counts'])

# convert the column names to those previously extracted
data.columns = names

# remove the first row (column names)
data = data.iloc[1:]

# remove the last row
data.drop(index=data.index[-1], axis = 0, inplace = True)

# remove the first column and last two columns
data = data.iloc[:, 1:12]

# rename columns to more useful names
data = data.rename(columns={'Identified Proteins (6388)':'Protein',
                            'Molecular Weight':'MW',
                            'Accession Number':'AN',
                            '4552_BirA-V5_1':'Control_BirAV5_1',
                            '4552_BirA-V5_2':'Control_BirAV5_2',
                            '4553_GFPmyc_1':'Control_GFPmyc_1',
                            '4553_GFPmyc_2':'Control_GFPmyc_2',
                            '4552_mTMEM-V5_1':'Treat_TMEMV5_1',
                            '4552_mTMEM-V5_2':'Treat_TMEMV5_2',
                            '4553_mTMEMmyc_1':'Treat_TMEMmyc_1',
                            '4553_mTMEMmyc_2':'Treat_TMEMmyc_2'})

# extract genenames to re-populate the Proteins column with only genenames
data['Protein'] = data['Protein'].str.extract(r'GN=([^ ]+)', expand = False)

# extract kDas out of MW column and convert to numeric
data['MW'] = pd.to_numeric(data['MW'].str.extract(r'(\d+)\s+kDa', expand = False))

# %%

# re-formatting the input to harmonize with the algorithm's requirements
#count_cols = [c for c in data.columns if c not in ["Protein"]]

# pivot longer
#long_data = data.melt(
#    id_vars=["Protein", "MW", "AN"],
#    value_vars=count_cols,
#    var_name="col",
#    value_name="count"
#)

# parse the condition/bait/rep column values
#def parse_col(col):
#    # Example: "CTRL_Bait1_rep2"
#    m = re.match(r"([^_]+)_([^_]+)_(\d+)", col)
#    if m:
#        condition, bait, rep = m.groups()
#        return condition, bait, int(rep)
#    return None, None, None

#long_data[["condition", "bait", "replicate"]] = long_data["col"].apply(
#    lambda x: pd.Series(parse_col(x))
#    )

#long_data = long_data.drop(columns=["col"])

#column_map = {
#    "bait": "bait",
#    "prey": "Protein",
#    "count": "count",
#    "condition": "condition",
#    "replicate": "replicate"
#}

#data.loc[data['Protein'].str.contains("GFP$", regex = True)]
#data.loc[data['Protein'].str.contains("TMEM184B", regex = True)]
# %% pipeline imports and metadata construction

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline
from saint.pipeline.classical_saint import run_classical_pipeline

metadata = {
    "biological_bait_names": {
        "TMEMV5": "TMEM184B",
        "TMEMmyc": "TMEM184B",
        "BirAV5": "birA",
        "GFPmyc": "GFP"
    },
    "AN" : {
        "TMEM184B" : "Q9Y519"
        },
    "MW" : {
        "TMEM184B" : 46
        }
}

# %% Run classical SAINT

results_classical = run_classical_pipeline(
    input_data=data,
    bait_names=['TMEMV5', 'TMEMmyc'],
    metadata=metadata,
    make_plots=True
    )

# %% Run hierarchical SAINT
results_hierarchical = run_hierarchical_pipeline(
    input_data=data,
    bait_names=["TMEMV5", "TMEMmyc"],
    metadata=metadata,
    make_plots=True
)

# %%

from saint.pipeline.hierarchical_saint import run_hierarchical_pipeline
from saint.diagnostics.diagnostics_tau_grid import diagnostics_tau_grid
from saint.diagnostics.diagnostics_hierarchical import make_hierarchical_plots
from saint.diagnostics.diagnostics_pipeline import plot_gamma3_density

em_results = results_hierarchical["raw_outputs"]["em_results"]
tau_info   = results_hierarchical["raw_outputs"]["tau_info"]

# per-bait hierarchical diagnostics
for bait, result in em_results.items():
    figs_hier = make_hierarchical_plots(result, bait)
    for f in figs_hier.values():
        f.show()

# per-bait tau-grid diagnostics
for bait, info in tau_info.items():
    diag_tau = diagnostics_tau_grid(info, bait)
    diag_tau["figure"].show()

# pipeline-level KDE
fig_kde = plot_gamma3_density(em_results)
fig_kde.show()



# %%

#output_df.loc[output_df['Protein'] == 'TMEM184B', ['Protein', 'lambda1', 'lambda2', 'lambda3', 'pi1', 'pi2', 'pi3', 'gamma1', 'gamma2', 'gamma3']]

# %% model validations

results_hierarchical["results_df"].head(20)
results_hierarchical["results_df"].sort_values("gamma3", ascending=False).head(20)

# %% Combined gamma3 density plot with TMEM184B marker

import seaborn as sns
import matplotlib.pyplot as plt

# Extract gamma3 values for each bait
df = results_hierarchical["results_df"]

gamma3_v5 = df[df["Bait"] == "TMEMV5"][["Protein", "gamma3"]].rename(columns={"gamma3": "gamma3_v5"})
gamma3_myc = df[df["Bait"] == "TMEMmyc"][["Protein", "gamma3"]].rename(columns={"gamma3": "gamma3_myc"})

# Merge and compute the average gamma3
merged = gamma3_v5.merge(gamma3_myc, on="Protein", how="inner")
merged["gamma3_avg"] = (merged["gamma3_v5"] + merged["gamma3_myc"]) / 2.0

# Extract TMEM184B average gamma3
tmem184b_value = float(merged.loc[merged["Protein"] == "TMEM184B", "gamma3_avg"])

# Plot
sns.set_style("whitegrid")
fig = plt.figure(figsize=(8, 5))

sns.kdeplot(merged["gamma3_avg"], fill=True, color="steelblue")
plt.axvline(tmem184b_value, color="red", linewidth=2)

plt.xlabel("Average gamma3 value")
plt.ylabel("Density")
plt.title("Combined gamma3 density for TMEMV5 and TMEMmyc with TMEM184B marker")

plt.show()

import os

fig.savefig("C:/Users/Erik/Desktop/Programming/R/Bio/Endo/analysis/modeling/pre_tuned_gamma3_density.png", dpi=300, bbox_inches="tight")
fig.savefig("C:/Users/Erik/Desktop/Programming/R/Bio/Endo/analysis/modeling/pre_tuned_gamma3_density.pdf", bbox_inches="tight")


# %% Ranked interactome table (cleaned, sorted, gamma3-based, no lambdas)

import pandas as pd

df = results["results_df"].copy()

# Remove unwanted columns if present
cols_to_drop = [
    c for c in ["pi", "tau", "gamma1", "gamma2", "lambda1", "lambda2", "lambda3"]
    if c in df.columns
]
df = df.drop(columns=cols_to_drop)

# Merge in gamma3_avg and FDR from the unique-protein table (merged)
df = df.merge(
    merged[["Protein", "gamma3_avg", "FDR"]],
    on="Protein",
    how="left"
)

# Rank proteins by gamma3_avg (1 = strongest interactor)
df["rank"] = df["gamma3_avg"].rank(method="dense", ascending=False).astype(int)

# Sort by rank, then by protein, then by bait
df = df.sort_values(["rank", "Protein", "Bait"])

# Reorder columns for clarity
ordered_cols = (
    ["rank", "Protein", "Bait"] +
    [c for c in df.columns if c not in ["rank", "Protein", "Bait"]]
)
df = df[ordered_cols]

ranked_interactome = df.reset_index(drop=True)

# If you want to save:
# ranked_interactome.to_csv("ranked_interactome.csv", index=False)

# %% Candidate interactors

# import matplotlib.pyplot as plt

# df = results["results_df"]

# gamma3_v5 = df[df["Bait"] == "TMEMV5"][["Protein", "gamma3"]].rename(columns={"gamma3": "gamma3_v5"})
# gamma3_myc = df[df["Bait"] == "TMEMmyc"][["Protein", "gamma3"]].rename(columns={"gamma3": "gamma3_myc"})

# merged = gamma3_v5.merge(gamma3_myc, on="Protein", how="inner")
# merged["gamma3_avg"] = (merged["gamma3_v5"] + merged["gamma3_myc"]) / 2.0

# cutoff = 0.95

# hits = merged[merged["gamma3_avg"] >= cutoff].sort_values("gamma3_avg", ascending=False)

# fig = plt.figure(figsize=(10, 5))
# plt.bar(hits["Protein"], hits["gamma3_avg"])
# plt.xticks(rotation=90)
# plt.ylabel("Average gamma3")
# plt.title(f"Proteins with average gamma3 at least {cutoff}")
# plt.tight_layout()
# plt.show()

# %% Compute principled FDR-based cutoff for gamma3_avg

import numpy as np
import pandas as pd

df = results["results_df"]

gamma3_v5 = df[df["Bait"] == "TMEMV5"][["Protein", "gamma3"]].rename(columns={"gamma3": "gamma3_v5"})
gamma3_myc = df[df["Bait"] == "TMEMmyc"][["Protein", "gamma3"]].rename(columns={"gamma3": "gamma3_myc"})

merged = gamma3_v5.merge(gamma3_myc, on="Protein", how="inner")
merged["gamma3_avg"] = (merged["gamma3_v5"] + merged["gamma3_myc"]) / 2.0

# Sort by gamma3_avg descending
merged = merged.sort_values("gamma3_avg", ascending=False).reset_index(drop=True)

# Compute FDR curve
merged["FDR"] = (1 - merged["gamma3_avg"]).cumsum() / (np.arange(len(merged)) + 1)

# Choose cutoff at FDR <= 0.05
alpha = 0.05
cutoff_gamma3 = merged.loc[merged["FDR"] <= alpha, "gamma3_avg"].min()

cutoff_gamma3

# %% Horizontal boxplots for proteins above the gamma3 cutoff

# import matplotlib.pyplot as plt
# import pandas as pd

# # Filter proteins above the principled cutoff
# hits = merged[merged["gamma3_avg"] >= cutoff_gamma3].copy()

# # Extract raw counts from the original results_df
# df_raw = results["results_df"][["Protein", "Bait", "rep1", "rep2"]]

# # Pivot to get 4 raw values per protein
# pivot = df_raw.pivot_table(index="Protein", columns="Bait", values=["rep1", "rep2"])

# # Flatten MultiIndex columns
# pivot.columns = [f"{rep}_{bait}" for rep, bait in pivot.columns]

# # Merge with hits
# hits_raw = hits.merge(pivot, on="Protein", how="left")

# # Order proteins by gamma3_avg
# hits_raw = hits_raw.sort_values("gamma3_avg", ascending=False)

# # Prepare data for horizontal boxplots
# data = []
# labels = []

# for _, row in hits_raw.iterrows():
#     values = []
#     for col in ["rep1_TMEMV5", "rep2_TMEMV5", "rep1_TMEMmyc", "rep2_TMEMmyc"]:
#         if col in row and not pd.isna(row[col]):
#             values.append(row[col])
#     data.append(values)
#     labels.append(row["Protein"])

# # Plot horizontal boxplots
# fig = plt.figure(figsize=(10, len(labels) * 0.4))
# plt.boxplot(data, labels=labels, vert=False)
# plt.xlabel("Raw counts")
# plt.title("Raw TMEMtag counts for proteins above gamma3 cutoff")
# plt.tight_layout()
# plt.show()

# %% Heatmap of raw counts for proteins above gamma3 cutoff (corrected)

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Filter proteins above the principled cutoff
hits = merged[merged["gamma3_avg"] >= cutoff_gamma3].copy()

# Ensure TMEM184B is included if it exists
if "TMEM184B" in merged["Protein"].dropna().values:
    tmem_row = merged[merged["Protein"] == "TMEM184B"]
    hits = pd.concat([tmem_row, hits]).drop_duplicates(subset=["Protein"])

# Extract raw counts
df_raw = results["results_df"][["Protein", "Bait", "rep1", "rep2"]]

# Pivot to get 4 raw values per protein
pivot = df_raw.pivot_table(index="Protein", columns="Bait", values=["rep1", "rep2"])
pivot.columns = [f"{rep}_{bait}" for rep, bait in pivot.columns]

# Merge with hits
hits_raw = hits.merge(pivot, on="Protein", how="left")

# Order proteins by gamma3_avg (highest at top)
hits_raw = hits_raw.sort_values("gamma3_avg", ascending=False)

# Select only the raw count columns
heatmap_data = hits_raw[
    ["rep1_TMEMV5", "rep2_TMEMV5", "rep1_TMEMmyc", "rep2_TMEMmyc"]
]

# Set index to protein names
heatmap_data.index = hits_raw["Protein"]

# Plot heatmap
sns.set_style("whitegrid")
fig = plt.figure(figsize=(8, len(heatmap_data) * 0.4))

ax = sns.heatmap(
    heatmap_data,
    cmap="viridis",
    annot=True,
    fmt=".0f",
    cbar_kws={"label": "Raw counts"}
)

plt.title("Raw TMEMtag counts for proteins above gamma3 cutoff")
plt.xlabel("Replicate")
plt.ylabel("Protein")
plt.tight_layout()
plt.show()

# %% Volcano data preparation (deduplicated lambdas, flipped SNR)

import numpy as np
import pandas as pd

df = results["results_df"]

# Deduplicate lambdas: one row per protein
lambda_df = (
    df[["Protein", "lambda1", "lambda2", "lambda3"]]
    .drop_duplicates(subset=["Protein"])
    .reset_index(drop=True)
)

# Merge with gamma3_avg table (merged)
volcano = merged.merge(lambda_df, on="Protein", how="left")

# Compute inverted SNR: log2((lambda1 + lambda2) / lambda3)
volcano["SNR"] = np.log2(
    (volcano["lambda1"] + volcano["lambda2"] + 1e-8) /
    (volcano["lambda3"] + 1e-8)
)

# Identify top hits by gamma3_avg
top_hits = (
    volcano.sort_values("gamma3_avg", ascending=False)
    .head(10)["Protein"]
    .tolist()
)

# Ensure TMEM184B is included
if "TMEM184B" in volcano["Protein"].dropna().values:
    if "TMEM184B" not in top_hits:
        top_hits.append("TMEM184B")
        
# %% Clean volcano (no histograms, no contours)

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("whitegrid")
fig = plt.figure(figsize=(8, 6))

# Background points (semi-transparent)
sns.scatterplot(
    data=volcano,
    x="SNR",
    y="gamma3_avg",
    color="gray",
    alpha=0.25,
    s=40
)

# Highlighted proteins (opaque)
highlight = volcano[volcano["Protein"].isin(top_hits)]
sns.scatterplot(
    data=highlight,
    x="SNR",
    y="gamma3_avg",
    hue="Protein",
    palette="tab10",
    s=90,
    edgecolor="black",
    linewidth=0.5
)

plt.axhline(cutoff_gamma3, color="red", linestyle="--", linewidth=1)

plt.xlabel("log2( (lambda1 + lambda2) / lambda3)\nSNR (Inverted Signal:Noise Ratio)")
plt.ylabel("Average gamma3")
plt.title("Pre-tuned Hierarchical EM SNR Volcano Plot")

ax = plt.gca()

# Round tick labels to 2 decimals
ticks = ax.get_xticks()
ax.set_xticklabels([round(t, 2) for t in ticks])

# Remove offset text
ax.get_xaxis().get_offset_text().set_visible(False)


plt.tight_layout()
plt.show()


fig.savefig("C:/Users/Erik/Desktop/Programming/R/Bio/Endo/analysis/modeling/pre_tuned_SNR_volcano.png", dpi=300, bbox_inches="tight")
fig.savefig("C:/Users/Erik/Desktop/Programming/R/Bio/Endo/analysis/modeling/pre_tuned_SNR_volcano.pdf", bbox_inches="tight")

# %% Volcano with marginal histograms

# sns.set_style("whitegrid")

# g = sns.JointGrid(
#     data=volcano,
#     x="SNR",
#     y="gamma3_avg",
#     height=8
# )

# # Background points
# g.plot_joint(
#     sns.scatterplot,
#     color="gray",
#     alpha=0.25,
#     s=40
# )

# # Highlighted proteins
# highlight = volcano[volcano["Protein"].isin(top_hits)]
# sns.scatterplot(
#     data=highlight,
#     x="SNR",
#     y="gamma3_avg",
#     hue="Protein",
#     palette="tab10",
#     s=90,
#     edgecolor="black",
#     linewidth=0.5,
#     ax=g.ax_joint
# )

# # Marginals
# sns.histplot(volcano["SNR"], ax=g.ax_marg_x, color="gray", alpha=0.4)
# sns.histplot(volcano["gamma3_avg"], ax=g.ax_marg_y, color="gray", alpha=0.4, orientation="horizontal")

# g.ax_joint.set_xlabel("log2( (lambda1 + lambda2) / lambda3)\nSNR (Inverted Signal:Noise Ratio)")
# g.ax_joint.set_ylabel("Average gamma3")
# g.ax_joint.set_title("Volcano Plot with Marginal Histograms")

# plt.tight_layout()
# plt.show()

# %% Volcano with density contours

# sns.set_style("whitegrid")
# fig = plt.figure(figsize=(8, 6))

# # KDE contours
# sns.kdeplot(
#     data=volcano,
#     x="SNR",
#     y="gamma3_avg",
#     fill=False,
#     levels=8,
#     color="gray",
#     alpha=0.5
# )

# # Background points
# sns.scatterplot(
#     data=volcano,
#     x="SNR",
#     y="gamma3_avg",
#     color="gray",
#     alpha=0.25,
#     s=40
# )

# # Highlighted proteins
# highlight = volcano[volcano["Protein"].isin(top_hits)]
# sns.scatterplot(
#     data=highlight,
#     x="SNR",
#     y="gamma3_avg",
#     hue="Protein",
#     palette="tab10",
#     s=90,
#     edgecolor="black",
#     linewidth=0.5
# )

# plt.axhline(cutoff_gamma3, color="red", linestyle="--", linewidth=1)

# plt.xlabel("log2( (lambda1 + lambda2) / lambda3)\nSNR (Inverted Signal:Noise Ratio)")
# plt.ylabel("Average gamma3")
# plt.title("Volcano Plot with Density Contours")
# plt.tight_layout()
# plt.show()

# %% Combined Volcano + Staggered Heatmap (TMEM + top 25 + bottom 10, raw counts)

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sns.set_style("whitegrid")

# -----------------------------
# Volcano prep
# -----------------------------
df = results["results_df"]

lambda_df = (
    df[["Protein", "lambda1", "lambda2", "lambda3"]]
    .drop_duplicates(subset=["Protein"])
    .reset_index(drop=True)
)

volcano = merged.merge(lambda_df, on="Protein", how="left")

volcano["SNR"] = np.log2(
    (volcano["lambda1"] + volcano["lambda2"] + 1e-8) /
    (volcano["lambda3"] + 1e-8)
)

top_hits = (
    volcano.sort_values("gamma3_avg", ascending=False)
    .head(10)["Protein"]
    .tolist()
)

if "TMEM184B" in volcano["Protein"].dropna().values:
    if "TMEM184B" not in top_hits:
        top_hits.append("TMEM184B")

highlight = volcano[volcano["Protein"].isin(top_hits)]

# -----------------------------
# Heatmap prep (TMEM + top 25 + bottom 10)
# -----------------------------
df_raw = df[["Protein", "Bait", "rep1", "rep2"]]

pivot = df_raw.pivot_table(index="Protein", columns="Bait", values=["rep1", "rep2"])
pivot.columns = [f"{rep}_{bait}" for rep, bait in pivot.columns]

hits = merged.copy()

tmem_rows = df[df["Protein"].str.contains("TMEM184B", case=False, na=False)][["Protein"]].drop_duplicates()
hits = pd.concat([tmem_rows, hits]).drop_duplicates(subset=["Protein"])

hits_raw = hits.merge(pivot, on="Protein", how="left")

# Identify TMEM block
tmem_block = hits_raw[hits_raw["Protein"].str.contains("TMEM184B", case=False, na=False)]

# Identify top 25 block
non_tmem = hits_raw[~hits_raw["Protein"].isin(tmem_block["Protein"])]
top25_block = non_tmem.sort_values("gamma3_avg", ascending=False).head(25)

# Identify bottom 10 block
remaining = non_tmem[~non_tmem["Protein"].isin(top25_block["Protein"])]
bottom10_block = remaining.sort_values("gamma3_avg", ascending=True).head(10)

# Extract raw counts
cols = ["rep1_TMEMV5", "rep2_TMEMV5", "rep1_TMEMmyc", "rep2_TMEMmyc"]

tmem_data = tmem_block[cols]
top25_data = top25_block[cols]
bottom10_data = bottom10_block[cols]

# Insert NaN rows as visual breaks
gap = pd.DataFrame([[np.nan]*len(cols)], columns=cols)

staggered = pd.concat([
    tmem_data,
    gap,
    top25_data,
    gap,
    bottom10_data
])

# Build y‑tick labels
labels = (
    list(tmem_block["Protein"]) +
    [""] +
    list(top25_block["Protein"]) +
    [""] +
    list(bottom10_block["Protein"])
)

# -----------------------------
# Combined figure
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 12), gridspec_kw={"width_ratios": [1.2, 1]})

# ---- Volcano (left) ----
ax = axes[0]

sns.scatterplot(
    data=volcano,
    x="SNR",
    y="gamma3_avg",
    color="gray",
    alpha=0.25,
    s=40,
    ax=ax
)

sns.scatterplot(
    data=highlight,
    x="SNR",
    y="gamma3_avg",
    hue="Protein",
    palette="tab10",
    s=90,
    edgecolor="black",
    linewidth=0.5,
    ax=ax
)

ax.axhline(cutoff_gamma3, color="red", linestyle="--", linewidth=1)
ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

ax.set_xlabel("log2( (lambda1 + lambda2) / lambda3)\nSNR (Inverted Signal:Noise Ratio)")
ax.set_ylabel("Average gamma3")
ax.set_title("Pre-tuned Hierarchical EM SNR Volcano Plot")



# ---- Staggered Heatmap (right) ----
ax2 = axes[1]

hm = sns.heatmap(
    staggered,
    cmap="viridis",
    annot=True,          # <-- raw values inside cells
    fmt=".0f",           # <-- integer formatting
    cbar_kws={"orientation": "vertical"},
    ax=ax2
)

# Horizontal title above colorbar
cbar = hm.collections[0].colorbar
cbar.ax.set_title("Raw counts", fontsize=10, pad=10)

# Explicit y‑ticks
ax2.set_yticks(np.arange(len(labels)) + 0.5)
ax2.set_yticklabels(labels, fontsize=9)   # <-- larger y-axis text

ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=9)

ax2.set_title("Raw TMEM-tag Counts\n(Staggered: TMEM → Top 25 → Bottom 10)")
ax2.set_xlabel("Replicate")
ax2.set_ylabel("Protein")

plt.tight_layout()
plt.show()

fig.savefig("C:/Users/Erik/Desktop/Programming/R/Bio/Endo/analysis/modeling/pre_tuned_volcano_heatmap.png", dpi=300, bbox_inches="tight")
fig.savefig("C:/Users/Erik/Desktop/Programming/R/Bio/Endo/analysis/modeling/pre_tuned_volcano_heatmap.pdf", bbox_inches="tight")