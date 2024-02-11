log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
spots.pass <- snakemake@input[["spots_pass"]]
seurat.object <- snakemake@input[["seurat"]]

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Slide
slide <- str_remove_all(str_extract(basename(snakemake@output[["raw_counts"]]),
                                    pattern = "_slide[0-9]{1}_"), pattern = "_")

# Spots after filtering
spots <- scan(spots.pass, sep = "\n", what = character())

# Spots in this slide
spots.slide.pass <- grep(pattern = slide, x = spots, value = TRUE)

# Raw counts
counts <- as.matrix(seuratobj@assays$Spatial@counts)[, spots.slide.pass]

# Coordinates
coordinates <- seuratobj@images[[slide]]@coordinates[spots.slide.pass, ] %>%
  select(col, row)

# Save matrices
write.table(counts, file = snakemake@output[["raw_counts"]],
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(coordinates, file = snakemake@output[["coordinates"]],
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
            
# Session info
sessionInfo()