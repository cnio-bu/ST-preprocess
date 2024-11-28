import sys
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = open(snakemake.log[0], "a+")

import random
import os
import scanpy as sc
import pandas as pd
import session_info
from scipy.io import mmread

## SNAKEMAKE I/O ##
mtx_file = snakemake.params["matrix"]
bc_file = snakemake.params["barcodes"]
features_file = snakemake.params["features"]
meta_file = snakemake.params["metadata"]
coord_file = snakemake.params["coordinates"]
gene_column = snakemake.params["gene_column"][0] - 1

## CODE ##
# Random seed
random.seed(1)

# Slide and patient
slide = os.path.basename(os.path.dirname(mtx_file))
patient = os.path.basename(os.path.dirname(os.path.dirname(mtx_file))).split("_")[0]

# Read 10X matrix
adata = sc.read_mtx(mtx_file)
adata_bc = pd.read_csv(bc_file, header = None)
adata_features = pd.read_csv(features_file, header = None, sep = "\t")
adata = adata.T
adata.obs["barcode"] = [x + "_" + slide for x in adata_bc[0].tolist()]
adata.var["gene_name"] = adata_features[gene_column].tolist()
adata.var.index = adata.var["gene_name"]
adata.var_names_make_unique()
adata.var.drop(columns=["gene_name"], inplace = True)

# Metadata
adata_meta = pd.read_csv(meta_file, header = 0)
adata_meta.rename(columns = { adata_meta.columns[0]: "barcode" }, inplace = True)
adata_meta["barcode"] = adata_meta["barcode"] + "_" + slide
adata_meta["patient"] = patient
adata_meta["slide"] = slide
spots_metadata = adata_meta["barcode"].tolist()
adata = adata[adata.obs["barcode"].isin(spots_metadata), :]
adata.obs = adata.obs.merge(adata_meta, on = "barcode")
adata.obs.set_index("barcode", inplace = True)

# Coordinates
spatial_coords = pd.read_csv(coord_file, header = None)
spatial_coords.columns = ["barcode", "in_tissue", "array_row", "array_col", "x", "y"]
spatial_coords["barcode"] = spatial_coords["barcode"] + "_" + slide
spatial_coords = spatial_coords[spatial_coords["in_tissue"] == 1]
spatial_coords.set_index("barcode", inplace = True)
spatial_coords = spatial_coords.reindex(adata.obs.index)
adata.obsm["spatial"] = spatial_coords[["x", "y"]].values

# Save object
adata.write_h5ad(str(snakemake.output))

# Session info
session_info.show()