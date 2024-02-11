log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(estimate))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]
out.dir <- dirname(snakemake@output[[1]])

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Normalised data
norm.counts <- as.matrix(seuratobj@assays$SCT@data)
norm.tsv <- paste0(out.dir, "/SCT_counts.tsv")
write.table(norm.counts, file = norm.tsv, sep = "\t", row.names = TRUE, 
            quote = FALSE)

# Common genes
commongenes.gct <- paste0(out.dir, "/SCT_counts_filtered.gct")
filterCommonGenes(input.f = norm.tsv, id = "GeneSymbol", 
                  output.f = commongenes.gct)

# ESTIMATE scores
estimate.gct <- paste0(out.dir, "/Scores.gct")
estimateScore(commongenes.gct, estimate.gct, platform = "illumina")

# Reformat ESTIMATE output
estimate <- t(read.table(file = estimate.gct, skip = 2, header = TRUE)) %>% 
  as.data.frame
colnames(estimate) <- estimate[1, ]
estimate <- estimate[-c(1:2), ] %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  mutate(ESTIMATEScaled = scales::rescale(ESTIMATEScore, to = c(0, 1)))
rownames(estimate) <- gsub("\\.", "-", rownames(estimate))

# Save table
write.table(estimate, file = snakemake@output[[1]], sep = "\t",
            row.names = TRUE, quote = FALSE)
            
# Session info
sessionInfo()