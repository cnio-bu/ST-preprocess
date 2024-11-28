log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
card.deconvolution <- unlist(snakemake@input)

## CODE ##
# Random seed
set.seed(1)

# Read deconvolution results and append them
if (length(card.deconvolution) > 1) {
    deconvolution <- lapply(card.deconvolution, FUN = function(x) {
        read.table(x, header = TRUE, row.names = 1, sep = "\t")
    }) |>
    bind_rows()
} else {
    deconvolution <-  read.table(card.deconvolution, header = TRUE, row.names = 1, sep = "\t")
}

# Save merged table
write.table(deconvolution, file = snakemake@output[[1]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
            
# Session info
sessionInfo()