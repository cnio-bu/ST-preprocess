log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(spacexr))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
raw.file <- snakemake@input[["raw_counts"]]
meta.file <- snakemake@input[["metadata"]]

## SNAKEMAKE PARAMETERS ##
lvl <- snakemake@params[["level"]]

## CODE ##
# Random seed
set.seed(1)

# Read raw counts
rawcounts <- read.table(raw.file, header = TRUE, row.names = 1, sep = "\t")

# Read metadata
metadata <- read.table(meta.file, header = TRUE, row.names = 1, sep = "\t")

# Extract cell types
cell.types.major <- metadata %>%
  mutate(celltype_major = as.factor(celltype_major)) %>%
  pull(celltype_major) %>%
  setNames(rownames(metadata))

# Create a reference object
reference.major <- Reference(rawcounts, cell_types = cell.types.major, n_max_cells = 500)

# Save reference
saveRDS(reference.major, file = snakemake@output[[1]])

# Session info
sessionInfo()