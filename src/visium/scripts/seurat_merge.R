log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.objects <- sort(unlist(snakemake@input))

## CODE ##
# Random seed
set.seed(1)

# Read Seurat objects
seuratobjs <- lapply(unlist(seurat.objects), FUN = function(x) readRDS(x))

# If there is more than one slide...
if (length(seuratobjs) > 1) {
    # Get common genes between all objects
    common.genes <- Reduce(intersect, 
                           lapply(seuratobjs, FUN = function(x) rownames(x)))

    # Filter objects to keep only common genes
    seuratobjs <- lapply(seuratobjs, FUN = function(x) x[common.genes, ])

    # Merge Seurat objects
    seuratobj.merged <- merge(seuratobjs[[1]], 
                              y = seuratobjs[2:length(seurat.objects)])
} else {
    # If there is one slide, return original Seurat object
    seuratobj.merged <- seuratobjs[[1]]
}

# Rename spots
patient <- unique(seuratobj.merged@meta.data$patient)
spot.names <- paste(colnames(seuratobj.merged), patient, sep = "_")
seuratobj.merged <- RenameCells(seuratobj.merged, new.names = spot.names)

# Save object
saveRDS(seuratobj.merged, file = snakemake@output[[1]])

# Session info
sessionInfo()