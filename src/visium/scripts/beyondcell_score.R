log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(beyondcell))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[["seurat"]]
gs.gmt <- snakemake@input[["geneset"]]
out.file <- basename(snakemake@output[[1]])

## SNAKEMAKE PARAMETERS ##
expr.thres <- snakemake@params[["expr_thres"]]

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Set assay
DefaultAssay(seuratobj) <- "SCT"

# Get signature collection
gs <- GenerateGenesets(gs.gmt, perform.reversal = FALSE)

# Compute BCS
bcobj <- suppressWarnings(
  bcScore(seuratobj, gs = gs, expr.thres = expr.thres)
)

# Number of NAs
n.NA <- data.frame(nNAs = colSums(is.na(bcobj@normalized)),
                   row.names = colnames(bcobj@normalized))
bcobj <- bcAddMetadata(bcobj, n.NA)
bcobj.recomputed <- bcobj

# Save object
saveRDS(bcobj.recomputed, file = snakemake@output[[1]])

# Session info
sessionInfo()