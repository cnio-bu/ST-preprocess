log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(CARD))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
reference <- snakemake@input[["ref"]]
raw.counts <- snakemake@input[["raw_counts"]]
coordinates <- snakemake@input[["coordinates"]]

## CODE ##
# Random seed
set.seed(1)

# Read spatial raw counts and coordinates
st.rawcounts <- read.table(raw.counts, header = TRUE, row.names = 1, sep = "\t") |>
  as.matrix()
colnames(st.rawcounts) <- 
  str_replace(colnames(st.rawcounts), pattern = "\\.", replacement = "-")

coords <- read.table(coordinates, header = TRUE, row.names = 1, sep = "\t") %>%
  rename(x = col, y = row)

# Read reference object
reference <- readRDS(reference)

# Get single-cell raw counts and metadata
sc.rawcounts <- as.matrix(reference@assays$RNA@counts)
sc.metadata <- reference@meta.data

# Create a CARD object
CARDobj <- 
  createCARDObject(sc_count = sc.rawcounts, sc_meta = sc.metadata, 
                   spatial_count = st.rawcounts, spatial_location = coords,
                   ct.varname = "celltype_major", 
                   ct.select = unique(sc.metadata$celltype_major), 
                   sample.varname = "Patient", minCountGene = 100, 
                   minCountSpot = 5)

# Spot deconvolution
CARDobj <- CARD_deconvolution(CARD_object = CARDobj)

# Cell type proportions for each spot
cell.prop <- CARDobj@Proportion_CARD

# Save cell proportions
write.table(cell.prop, file = snakemake@output[[1]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Session info
sessionInfo()