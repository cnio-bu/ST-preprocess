log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(SCEVAN))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]
out.dir <- dirname(snakemake@output[[1]])

## CODE ##
# Random seed
set.seed(1)

# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Raw data
raw.counts <- as.matrix(seuratobj@assays$Spatial@data)

# Reference
reference <- seuratobj@meta.data %>%
  filter(Tumour == "Non-tumour") %>%
  rownames()

# Run SCEVAN
setwd(out.dir)
results <- pipelineCNA(raw.counts, norm_cell = reference, beta_vega = 0.5,
                       SUBCLONES = TRUE, ClonalCN = TRUE, plotTree = TRUE, 
                       SCEVANsignatures = TRUE, organism = "human")

# Get subclones
subclones <- results %>%
  select(subclone)

# Add metadata
seuratobj <- AddMetaData(seuratobj, subclones)

# Relabel spots detected as subclones as "Tumour"
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(subclone = case_when(is.na(subclone) & Tumour == "Tumour" ~ "diploid",
                              is.na(subclone) & Tumour != "Tumour" ~ NA_character_,
                              TRUE ~ as.character(subclone)),
         Tumour = if_else(!is.na(subclone), "Tumour", "Non-tumour"))

# Save object
saveRDS(seuratobj, file = snakemake@output[[1]])

# Session info
sessionInfo()