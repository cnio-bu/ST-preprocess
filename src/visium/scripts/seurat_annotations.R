log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]
out.dir <- dirname(snakemake@output[["seurat"]])

### FUNCTIONS ###
# Returns Seurat's SpatialFeaturePlot default colourscale
colourscale <- function(n = 100, rev = TRUE) {
  colours <- RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  if (rev) {
    colours <- colours |> 
      rev()
  }
  colours <- colours |>
    grDevices::colorRampPalette()
  return(colours(n = n))
}

# Draws custom SpatialFeaturePlots
customSpatialFeaturePlot <- function(x, features, limits = NULL, 
                                     legend.title = NULL, ...) {
  if (is.null(limits)) {
    scale <- legend <- titles <- NULL
  } else {
    scale <- scale_fill_gradientn(colours = colourscale(n = 100, ...),
                                  limits = limits)
    legend <- labs(fill = legend.title)
    titles <- rep(features, each = n.slides)
  }
  lplots <- SpatialFeaturePlot(x, features = features, pt.size.factor = 2, 
                               combine = FALSE, image.alpha = 0)
  lplots <- suppressMessages(
    lapply(seq_along(lplots), FUN = function(i) {
      lplots[[i]] + 
        scale + 
        ggtitle(titles[[i]]) +
        legend
    }))
  return(lplots)
}

# Draws custom SpatialDimPlots
customSpatialDimPlot <- function(x, group) {
  lplots <- SpatialDimPlot(x, group.by = group, pt.size.factor = 2, 
                           combine = FALSE, image.alpha = 0)
  return(lplots)
}

## CODE ##
# Random seed
set.seed(1)

# Create outdirs
dir.create(paste0(out.dir, "/plots/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(out.dir, "/ggplots/"), showWarnings = FALSE)

# Rename outdir
file.rename(paste0(out.dir, "/output/"), paste0(out.dir, "/clones/"))

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Plot pathologist annotations
slides <- seuratobj@images
n.slides <- length(slides)
is.empty.patho <- all(is.na(seuratobj@meta.data$Pathologist))
if (is.empty.patho) {
  pathologist <- lapply(seq_along(slides), FUN = function(i) plot_spacer())
} else {
  pathologist <- 
    customSpatialDimPlot(seuratobj, group = "Pathologist")
  ggsave(wrap_plots(pathologist, nrow = 1), width = 9 * n.slides, 
         filename = paste0(out.dir, "/plots/pathologist_annotations.pdf"))
}

# Plot cell type proportions
cell.types <- seuratobj@meta.data %>%
  select(B.cells:T.cells) %>%
  colnames()
deconvolution.list <- 
  customSpatialFeaturePlot(seuratobj, limits = c(0, 1),
                           features = cell.types,
                           legend.title = "Cell Proportion")
                          
deconvolution.cancer <- 
  customSpatialFeaturePlot(seuratobj, limits = c(0, 1),
                           features = "Cancer.Epithelial",
                           legend.title = "Cell Proportion")

# Plot ESTIMATE annotations
estimate.stromal <- 
  customSpatialFeaturePlot(seuratobj, feature = "StromalScore")
estimate.immune <- 
  customSpatialFeaturePlot(seuratobj, feature = "ImmuneScore")
estimate.purity <- 
  customSpatialFeaturePlot(seuratobj, feature = "ESTIMATEScaled")
estimate.purity <- lapply(estimate.purity, FUN = function(x) {
  x + scale_fill_gradientn(colours = colourscale(n = 100, rev = FALSE), 
                           limits = c(0, 1))
})

# Plot aggregated annotation
tumour.tme <- customSpatialDimPlot(seuratobj, group = "Tumour")

# Patchworks
for (i in seq_along(slides)) {
  idx <- seq(i, length(deconvolution.list), by = n.slides)
  deconvolution <- wrap_plots(deconvolution.list[idx], nrow = 3, ncol = 3, 
                              guides = "collect") & theme(legend.position = "bottom")
  estimate <- estimate.stromal[[i]] | estimate.immune[[i]] | estimate.purity[[i]]
  complete.anno <- (pathologist[[i]] | estimate.purity[[i]]) / 
    (deconvolution.cancer[[i]] | tumour.tme[[i]])
  ggsave(deconvolution, width = 25, height = 35,
         filename = paste0(out.dir, "/plots/cell_proportions_slide", i, ".pdf"))
  ggsave(estimate, width = 25,
         filename = paste0(out.dir, "/plots/estimate_annotations_slide", i, ".pdf"))
  ggsave(complete.anno, width = 25, height = 20,
         filename = paste0(out.dir, "/plots/tumour_tme_annotations_slide", i, ".pdf"))
}

# Plot subclones
subclones <- customSpatialDimPlot(seuratobj, group = "subclone")
ggsave(wrap_plots(subclones), width = 8 * n.slides,
         filename = paste0(out.dir, "/plots/subclones.pdf"))

# Save plots
all.plots <- list(pathologist, deconvolution.list, estimate.stromal, 
                  estimate.immune, estimate.purity, tumour.tme, subclones)
save(all.plots, file = paste0(out.dir, "/ggplots/seurat_annotations.RData"))

# Transposed SCT counts
counts <- t(as.matrix(seuratobj@assays$SCT@data))

# Metadata
metadata <- seuratobj@meta.data

# Save matrices
write.table(counts, file = snakemake@output[["sct_counts"]],
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(metadata, file = snakemake@output[["metadata"]],
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

# Save object
saveRDS(seuratobj, file = snakemake@output[["seurat"]])

# Session info
sessionInfo()