log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
mtx.dir <- dirname(snakemake@input[["matrix"]])
meta.file <- snakemake@input[["metadata"]]

## CODE ##
# Random seed
set.seed(1)

# Read 10X matrix
mtx <- Read10X(data.dir = mtx.dir)

# Read metadata
meta <- read.csv(meta.file, row.names = 1)
meta <- meta[-1, ] %>%
  mutate(nCount_RNA = as.numeric(nCount_RNA),
         nFeature_RNA = as.numeric(nFeature_RNA),
         Percent_mito = as.numeric(Percent_mito))

# Subset matrix and metadata
intersect.cells <- intersect(colnames(mtx), rownames(meta))
mtx <- mtx[, intersect.cells]
meta <- meta[intersect.cells, ]

# Create Seurat object
seuratobj <- CreateSeuratObject(counts = mtx, assay = "RNA", 
                                meta.data = meta, project = "Single-cell")

# Save object
saveRDS(seuratobj, file = snakemake@output[[1]])

# Session info
sessionInfo()