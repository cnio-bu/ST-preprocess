log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(Giotto))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]

## SNAKEMAKE PARAMETERS ##
regress.cell.cycle <- snakemake@params[["regress_cell_cycle"]]
idents <- snakemake@params[["idents"]]
lvl <- snakemake@params[["level"]]

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Split Seurat object
seuratobjs <- SplitObject(seuratobj, split.by = idents)

# Select desired level
seuratobj.splitted <- seuratobjs[[lvl]]

# Convert to Giotto object
giottoobj <- seuratToGiottoV4(seuratobj.splitted, spatial_assay = "Spatial")

# Normalise counts using the standard method, z-scoring feats over cells
giottoobj <- 
  normalizeGiotto(gobject = giottoobj, scalefactor = 6000, verbose = TRUE)

# Cell cycle effect
normcounts <- 
  getExpression(giottoobj, values = "normalized", output = "matrix", 
                set_defaults = TRUE) |>
  as.matrix()

g2m.genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% rownames(normcounts)]
s.genes <- cc.genes$s.genes[cc.genes$s.genes %in% rownames(normcounts)]
g2m.score <- colMeans(normcounts[g2m.genes, , drop = FALSE], na.rm = TRUE)
s.score <- colMeans(normcounts[s.genes, , drop = FALSE], na.rm = TRUE)

meta <- data.frame(G2M.Score = g2m.score, S.Score = s.score) %>%
  rownames_to_column("spot")

giottoobj <- 
  addCellMetadata(giottoobj, spat_unit = "cell", feat_type = "rna", 
                  new_metadata = meta, by_column = TRUE, column_cell_ID = "spot")

# Regress region and cell cycle, if needed
vars.regress <- c()

if (regress.cell.cycle) {
    vars.regress <- c(vars.regress, c("S.Score", "G2M.Score"))
}

n.regions <- seuratobj.splitted@meta.data %>%
  pull(region) %>%
  unique() %>%
  length()

if (n.regions > 1) {
    region <- pDataDT(giottoobj)$region
    region.dummies <- model.matrix(~region - 1) %>%
      as.data.frame() %>%
      mutate(spot = pDataDT(giottoobj)$cell_ID)
    giottoobj <- 
      addCellMetadata(giottoobj, spat_unit = "cell", feat_type = "rna", 
                      new_metadata = region.dummies, by_column = TRUE, column_cell_ID = "spot")
    vars.regress <- c(vars.regress, colnames(region.dummies)[colnames(region.dummies) != "spot"])
}

if (length(vars.regress) > 0) {
  giottoobj <- 
    adjustGiottoMatrix(gobject = giottoobj, expression_values = c("normalized"),
                       covariate_columns = vars.regress)
  normcounts <- 
    getExpression(giottoobj, values = "normalized", output = "matrix", 
                  set_defaults = TRUE) |>
    as.matrix()
}

# Save Giotto object and normalised matrix
saveRDS(giottoobj, file = snakemake@output[["giotto"]])
write.table(normcounts, file = snakemake@output[["norm_counts"]],
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

# Session info
sessionInfo()