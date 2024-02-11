log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Clean pathologist annotation
pathologist <- seuratobj@meta.data %>%
  mutate(Pathologist = case_when(Classification %in% c("Uncertain", "") ~ 
                                   NA_character_,
                                 TRUE ~ Classification)) %>%
  select(Pathologist)

# Save table
write.table(pathologist, file = snakemake@output[[1]], sep = "\t",
            row.names = TRUE, quote = FALSE)
            
# Session info
sessionInfo()