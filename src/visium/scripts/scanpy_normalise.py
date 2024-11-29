import sys
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = open(snakemake.log[0], "a+")

import random
import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import session_info

## SNAKEMAKE I/O ##
cellcycle_genes = pd.read_csv(snakemake.input["cell_cycle"], sep = "\t")
scanpy_objects = list(snakemake.input["scanpy"])
with open(snakemake.input["spots_pass"]) as file:
    spots_pass = file.read().splitlines()

## SNAKEMAKE PARAMETERS ##
regress_cellcycle = snakemake.params["regress_cell_cycle"]

## CODE ##
# Random seed
random.seed(1)

# Read Scanpy objects
adatas = {}

for filename in scanpy_objects:
    slide = os.path.basename(filename).split("_")[1]
    adata = sc.read_h5ad(filename)
    adatas[slide] = adata

# Concatenate Scanpy objects
adata = ad.concat(adatas, label = "slide")

# Rename spots
patient = adata.obs["patient"].tolist()
patient = list(set(patient))[0]
adata.obs_names = adata.obs_names + "_" + patient

# Keep only pass spots
adata = adata[spots_pass].copy()

# Normalise, log-transform and scale the data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
normalized_counts_df = pd.DataFrame(
    adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,  # Convert sparse to dense if needed
    index = adata.obs_names,  # Cell names
    columns = adata.var_names  # Gene names
)
sc.pp.scale(adata)

# Cell cycle effect
s_genes = cellcycle_genes[cellcycle_genes["phase"] == "S"]["gene"].tolist()
g2m_genes = cellcycle_genes[cellcycle_genes["phase"] == "G2M"]["gene"].tolist()
sc.tl.score_genes_cell_cycle(adata, s_genes = s_genes, g2m_genes = g2m_genes)

# Regress region and cell cycle, if needed
vars_regress = []

if regress_cellcycle:
    vars_regress.extend(["S_score", "G2M_score"])

regions = set(adata.obs["region"].tolist())

if len(regions) > 1:
   region_dummies = pd.get_dummies(adata.obs["region"], prefix = "region")
   region_dummies = region_dummies.astype(int)
   adata.obs = pd.concat([adata.obs, region_dummies], axis = 1)
   vars_regress.extend(region_dummies.columns.tolist())
   
if len(vars_regress) > 0:
    sc.pp.regress_out(adata, vars_regress)
    adata.X = np.clip(adata.X, a_min = 0, a_max = None)  # Replace negative values with 0
    normalized_counts_df = pd.DataFrame(
        adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,  # Convert sparse to dense if needed
        index = adata.obs_names,  # Cell names
        columns = adata.var_names  # Gene names
    )

# Save normalised matrix
normalized_counts_df.T.to_csv(str(snakemake.output), sep = "\t")

# Session info
session_info.show()