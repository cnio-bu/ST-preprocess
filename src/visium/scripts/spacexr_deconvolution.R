log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(spacexr))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
reference <- snakemake@input[["ref"]]
raw.counts <- snakemake@input[["raw_counts"]]
coordinates <- snakemake@input[["coordinates"]]

## CODE ##
# Random seed
set.seed(1)

# Read raw counts and coordinates
rawcounts <- read.table(raw.counts, header = TRUE, row.names = 1, sep = "\t")
colnames(rawcounts) <- str_replace(colnames(rawcounts), pattern = "\\.",
                                   replacement = "-")

coords <- read.table(coordinates, header = TRUE, row.names = 1, sep = "\t")

# Read reference object
reference <- readRDS(reference)

# Create SpatialRNA object
puck <- SpatialRNA(coords, rawcounts)

# Run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 4, UMI_min = 100, 
                      CELL_MIN_INSTANCE = 25, MAX_MULTI_TYPES = 4)

myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

# Cell type proportions for each spot
cell.prop <- myRCTD@results$weights

# Normalise cell proportions to sum 1
cell.prop.normalised <- as.matrix(normalize_weights(cell.prop))

# Save cell proportions
write.table(cell.prop.normalised, file = snakemake@output[["deconvolution"]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
            
# Session info
sessionInfo()