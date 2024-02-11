log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]
out.dir <- dirname(snakemake@output[[1]])

## SNAKEMAKE PARAMETERS ##
regress.cell.cycle <- snakemake@params[["regress_cell_cycle"]]

## CODE ##
# Random seed
set.seed(1)

# Create patchwork theme
title_theme <- theme(plot.title = element_text(size = 16, face = "bold"))

# Create outdirs
dir.create(paste0(out.dir, "/cell_cycle/plots/without_regression"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(out.dir, "/cell_cycle/ggplots"), showWarnings = FALSE)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Regress region, if needed
vars.regress <- NULL
n.regions <- seuratobj@meta.data %>%
  pull(region) %>%
  unique() %>%
  length()
if (n.regions > 1) vars.regress <- "region"

# Perform SCT normalisation without cell cycle regression
seuratobj <- 
  SCTransform(seuratobj, assay = "Spatial", new.assay.name = "SCT",
              vars.to.regress = vars.regress,
              return.only.var.genes = TRUE, verbose = TRUE)

# Cell cycle effect
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
seuratobj.phase <- CellCycleScoring(seuratobj, g2m.features = g2m.genes, 
                                    s.features = s.genes, assay = "SCT")
Idents(seuratobj.phase) <- "Phase"

# Run PCA without cell cycle regression
seuratobj.phase <- RunPCA(seuratobj.phase, assay = "SCT", npcs = 50, 
                          features = VariableFeatures(seuratobj.phase))

# PCA plot without cell cycle regression
cell.cycle.without <- DimPlot(seuratobj.phase)
cell.cycle.phase.without <- DimPlot(seuratobj.phase, split.by = "Phase")

ggsave(cell.cycle.phase.without, width = 17,
       filename = paste0(out.dir, "/cell_cycle/plots/without_regression/", 
                         "cell_cycle_phase_without_regression.pdf"))

all.plots <- list(cell.cycle.without, cell.cycle.phase.without)

# If regress.cell.cycle is TRUE...
if (regress.cell.cycle) {
  # Create outdir
  dir.create(paste0(out.dir, "/cell_cycle/plots/with_regression"), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Regress by slide and cell cycle
  seuratobj.phase <- 
  SCTransform(seuratobj.phase, assay = "Spatial", new.assay.name = "SCT",
              return.only.var.genes = TRUE, verbose = TRUE,
              vars.to.regress = c(vars.regress, c("S.Score", "G2M.Score")))

  # Run PCA with cell cycle regression
  seuratobj.phase <- RunPCA(seuratobj.phase, assay = "SCT", npcs = 50, 
                            features = VariableFeatures(seuratobj.phase))
  
  # PCA plot with cell cycle regression
  cell.cycle.with <- DimPlot(seuratobj.phase)
  cell.cycle.phase.with <- DimPlot(seuratobj.phase, split.by = "Phase")
  
  ggsave(cell.cycle.phase.with, width = 17,
         filename = paste0(out.dir, "/cell_cycle/plots/with_regression/", 
                           "cell_cycle_phase_with_regression.pdf"))
  
  all.plots <- c(all.plots, list(cell.cycle.with, cell.cycle.phase.with))
  
  # Patchwork of cell cycle DimPlots
  cell.cycle <- (cell.cycle.without + ggtitle("Without regression") | 
                   cell.cycle.with + ggtitle("With regression")) + 
    plot_annotation(title = "Cell Cycle", theme = title_theme) + 
    plot_layout(guides = "collect")
  
  ggsave(cell.cycle, width = 14,
         filename = paste0(out.dir, "/cell_cycle/plots/cell_cycle.pdf"))
  
} else {
  # Save cell cycle DimPlot without cell cycle regression
  ggsave(cell.cycle.without, 
         filename = paste0(out.dir, "/cell_cycle/plots/without_regression/",
                           "cell_cycle_without_regression.pdf"))
}

# Save plots
save(all.plots, file = paste0(out.dir, "/cell_cycle/ggplots/cell_cycle.RData"))

# Save object
saveRDS(seuratobj.phase, file = snakemake@output[[1]])

# Session info
sessionInfo()