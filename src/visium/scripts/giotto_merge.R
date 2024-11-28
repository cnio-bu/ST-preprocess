log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
giotto.norm.counts  <- snakemake@input[["norm_counts"]]
giotto.deconvolution <- snakemake@input[["deconvolution"]]

## CODE ##
# Random seed
set.seed(1)

# Read normalised counts and append them
if (length(giotto.norm.counts) > 1) {
    normcounts <- lapply(giotto.norm.counts, FUN = function(x) {
        read.table(x, header = TRUE, row.names = 1, sep = "\t")
    }) |>
    bind_cols()
} else {
    normcounts <- read.table(giotto.norm.counts, header = TRUE, row.names = 1, sep = "\t")
}

# Read deconvolution results and append them
if (length(giotto.deconvolution) > 1) {
    deconvolution <- lapply(giotto.deconvolution, FUN = function(x) {
        read.table(x, header = TRUE, row.names = 1, sep = "\t")
    }) |>
    bind_rows()
} else {
    deconvolution <-  read.table(giotto.deconvolution, header = TRUE, row.names = 1, sep = "\t")
}

# Save merged tables
write.table(normcounts, file = snakemake@output[["norm_counts"]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(deconvolution, file = snakemake@output[["deconvolution"]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Session info
sessionInfo()