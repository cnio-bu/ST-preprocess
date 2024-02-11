log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[["seurat"]]
pathologist.file <- snakemake@input[["pathologist"]]
deconvolution.file <- snakemake@input[["deconvolution"]]
estimate.file <- snakemake@input[["estimate"]]
out.dir <- dirname(snakemake@output[[1]])

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Read pathologist metadata
pathologist.meta <- read.table(pathologist.file, header = TRUE, row.names = 1,
                               sep = "\t") %>%
  rownames_to_column("spot")

# Read deconvolution metadata
deconvolution.meta <- lapply(deconvolution.file, FUN = function(x) {
  read.table(x, header = TRUE, row.names = 1, sep = "\t") %>%
    rownames_to_column("spot")
}) %>%
  bind_rows()

# Read estimate metadata
estimate.meta <- read.table(estimate.file, header = TRUE, row.names = 1, 
                            sep = "\t") %>%
  rownames_to_column("spot")

# Format deconvolution metadata
deconvolution.long <- deconvolution.meta %>%
  pivot_longer(cols = colnames(deconvolution.meta)[-1], names_to = "cell", 
               values_to = "proportion") %>%
  mutate(cell = str_replace_all(cell, pattern = "\\.", replacement = " "),
         proportion = round(proportion, 2))

# Merge metadata
metadata <- pathologist.meta %>% 
  inner_join(deconvolution.meta) %>%
  inner_join(estimate.meta) %>%
  column_to_rownames("spot")

# Add metadata
seuratobj.anno <- AddMetaData(seuratobj, metadata)

# Create compartment annotation according to pathologist annotations and 
# ESTIMATE score
seuratobj.anno@meta.data <- seuratobj.anno@meta.data %>%
  mutate(Tumour.patho = Pathologist == "DCIS" | 
           str_detect(tolower(Classification), pattern = "cancer") |
           str_detect(tolower(Classification), pattern = "carcinoma"),
         Tumour.estimate = ESTIMATEScaled <= 0.4,
         Tumour.deconv = Cancer.Epithelial >= 0.6,
         Tumour = case_when(Tumour.patho & Tumour.estimate ~ "Tumour",
                            Tumour.patho & Tumour.deconv ~ "Tumour",
                            Tumour.estimate & Tumour.deconv ~ "Tumour",
                            TRUE ~ "Non-tumour")) %>%
  select(-c(Tumour.patho, Tumour.estimate, Tumour.deconv))

# Save object
saveRDS(seuratobj.anno, file = snakemake@output[[1]])

# Session info
sessionInfo()