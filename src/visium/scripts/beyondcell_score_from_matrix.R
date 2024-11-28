log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(beyondcell))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
scanpy.norm.counts <- snakemake@input[["scanpy_norm_counts"]]
giotto.norm.counts <- snakemake@input[["giotto_norm_counts"]]
gs.gmt <- snakemake@input[["geneset"]]

## SNAKEMAKE PARAMETERS ##
expr.thres <- snakemake@params[["expr_thres"]]

### FUNCTIONS ###
# Converts NA BCS to zero
NAto0 <- function(bc) {
  n.NA <- data.frame(nNAs = colSums(is.na(bc@normalized)),
                     row.names = colnames(bc@normalized))
  bcobj <- bcAddMetadata(bc, n.NA)
  return(bc)
}

## CODE ##
# Random seed
set.seed(1)

# Read normalised counts
scanpy.normcounts <- 
  read.table(scanpy.norm.counts, header = TRUE, row.names = 1, sep = "\t") %>%
  as.matrix()
colnames(scanpy.normcounts) <- 
  str_replace(colnames(scanpy.normcounts), pattern = "\\.", replacement = "-")
print(summary(as.vector(scanpy.normcounts)))

giotto.normcounts <- 
  read.table(giotto.norm.counts, header = TRUE, row.names = 1, sep = "\t") %>%
  as.matrix()
colnames(giotto.normcounts) <- 
  str_replace(colnames(giotto.normcounts), pattern = "\\.", replacement = "-")
print(summary(as.vector(giotto.normcounts)))

# Get signature collection
gs <- GenerateGenesets(gs.gmt, perform.reversal = FALSE)

# Compute BCS
scanpy.bcobj <- suppressWarnings(
  bcScore(scanpy.normcounts, gs = gs, expr.thres = expr.thres)
)

giotto.bcobj <- suppressWarnings(
  bcScore(giotto.normcounts, gs = gs, expr.thres = expr.thres)
)

# Convert NAs to 0
scanpy.bcobj.recomputed <- NAto0(scanpy.bcobj)
giotto.bcobj.recomputed <- NAto0(giotto.bcobj)

# Save object
saveRDS(scanpy.bcobj.recomputed, file = snakemake@output[["scanpy_bc"]])
saveRDS(giotto.bcobj.recomputed, file = snakemake@output[["giotto_bc"]])

# Session info
sessionInfo()