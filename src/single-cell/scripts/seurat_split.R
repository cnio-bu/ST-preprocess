log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]
out.dir <- dirname(snakemake@output[[1]])

## SNAKEMAKE PARAMETERS ##
idents <- snakemake@params[["idents"]]
lvl <- snakemake@params[["level"]]

## CODE ##
# Random seed
set.seed(1)

# Create outdirs
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Split Seurat object
seuratobjs <- SplitObject(seuratobj, split.by = idents)

# Select desired level
seuratobj.splitted <- seuratobjs[[lvl]]

# Raw counts
rawcounts <- GetAssayData(seuratobj.splitted, slot = "counts", assay = "RNA") |>
  as.matrix()

# Metadata
metadata <- seuratobj.splitted@meta.data %>%
  select(Patient, celltype_major:CNA_value)

# Save raw counts
write.table(rawcounts, file = snakemake@output[["raw_counts"]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Save metadata
write.table(metadata, file = snakemake@output[["metadata"]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Save object
saveRDS(seuratobj.splitted, file = snakemake@output[["seurat"]])

# Session info
sessionInfo()