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
mt.limits <- snakemake@params[["mt_limits"]]
rb.limits <- snakemake@params[["rb_limits"]]
count.limits <- snakemake@params[["count_limits"]]
feat.limits <- snakemake@params[["feat_limits"]]

### FUNCTIONS ###
# Computes min and max limits 
# Also converts to numeric any non min/max value
getMinMax <- function(limits, feature, obj) {
  is.min <- limits == "min"
  is.max <- limits == "max"
  if (any(is.min)) limits[is.min] <- min(obj@meta.data[[feature]])
  if (any(is.max)) limits[is.max] <- max(obj@meta.data[[feature]])
  if (any(c(is.min, is.max))) limits <- as.numeric(limits)
  return(limits)
}

# Checks that the threshold limits are in the right format and reorders them
# If both limits are NULL, returns NULL
# Requires to specify the metadata column (feature) and the seurat object
setLimits <- function(limits, feature, obj) {
  if (length(limits) > 0 | !is.null(limits)[1]) {
    ### Checks
    varname <- deparse(substitute(limits))
    if (length(limits) != 2) stop(paste(varname, 'must be of length 2.'))
    limits <- getMinMax(limits, feature, obj)
    if (!all(is.numeric(limits))) {
        stop(paste(varname, 'must be numeric or null.'))
    }
    if (any(limits < 0)) stop(paste(varname, 'must be >= 0.'))
    ### Code
    limits <- sort(limits)
  }
  return(limits)
}

# Draws custom VlnPlot
customVln <- function(x, features) {
  VlnPlot(x, features = features, pt.size = 0.1) + blank_x_theme
}

## CODE ##
# Random seed
set.seed(1)

# Create ggplot2 theme
blank_x_theme <- theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(), 
                       axis.ticks.x = element_blank(),
                       plot.title = element_text(size = 14, face = "plain"))

# Create patchwork theme
title_theme <- theme(plot.title = element_text(size = 16, face = "bold"))

# Create outdirs
dir.create(paste0(out.dir, "/plots"), recursive = TRUE, 
           showWarnings = FALSE)
dir.create(paste0(out.dir, "/ggplots"), showWarnings = FALSE)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Percentage ribo (mito is already computed)
seuratobj[["Percent_ribo"]] <- 
  PercentageFeatureSet(object = seuratobj, pattern = "^RP[SL][[:digit:]]")
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(Percent_ribo = if_else(is.na(Percent_ribo), 0, Percent_ribo))

# VlnPlots before cell filtering
subtype_cols <- scales::hue_pal()(3)
Idents(seuratobj) <- "Patient"
vln.mt.patient.before <- 
  VlnPlot(seuratobj, features = "Percent_mito", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Pre-filtering")
vln.rb.patient.before <- 
  VlnPlot(seuratobj, features = "Percent_ribo", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Pre-filtering")
vln.count.patient.before <- 
  VlnPlot(seuratobj, features = "nCount_RNA", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Pre-filtering")
vln.feat.patient.before <- 
  VlnPlot(seuratobj, features = "nFeature_RNA", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Pre-filtering")

Idents(seuratobj) <- "subtype"
vln.mt.tumour.before <- customVln(seuratobj, features = "Percent_mito") + 
  ggtitle("Pre-filtering")
vln.rb.tumour.before <- customVln(seuratobj, features = "Percent_ribo") +
  ggtitle("Pre-filtering")
vln.count.tumour.before <- customVln(seuratobj, features = "nCount_RNA") +
  ggtitle("Pre-filtering")
vln.feat.tumour.before <- customVln(seuratobj, features = "nFeature_RNA") +
  ggtitle("Pre-filtering")

Idents(seuratobj) <- "celltype_major"
vln.mt.cell.before <- customVln(seuratobj, features = "Percent_mito") + 
  ggtitle("Pre-filtering")
vln.rb.cell.before <- customVln(seuratobj, features = "Percent_ribo") +
  ggtitle("Pre-filtering")
vln.count.cell.before <- customVln(seuratobj, features = "nCount_RNA") +
  ggtitle("Pre-filtering")
vln.feat.cell.before <- customVln(seuratobj, features = "nFeature_RNA") +
  ggtitle("Pre-filtering")

# QC variables
qcvars <- c("Percent_mito", "Percent_ribo", "nCount_RNA", "nFeature_RNA")

# Set limits
limits.list <- list(mt.limits, rb.limits, count.limits, feat.limits)
limits.list <- lapply(seq_along(limits.list), FUN = function(i) {
  setLimits(limits.list[[i]], feature = qcvars[i], obj = seuratobj)
})
names(limits.list) <- qcvars

# Cell-level filtering
cells <- FetchData(seuratobj, vars = qcvars)
cells.subset <- lapply(names(limits.list), FUN = function(x) {
  if (is.null(limits.list[[x]])) {
    cells.subset <- cells %>% rownames()
  } else {
    cells.subset <- cells %>% 
      filter(get(x) >= limits.list[[x]][1] & get(x) <= limits.list[[x]][2]) %>%
      rownames()
  }
})
cells.pass <- Reduce(intersect, cells.subset)
seuratobj.filtered <- seuratobj[, cells.pass]

# VlnPlots after cell filtering
Idents(seuratobj.filtered) <- "Patient"
vln.mt.patient.after <- 
  VlnPlot(seuratobj.filtered, features = "Percent_mito", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Post-filtering")
vln.rb.patient.after <- 
  VlnPlot(seuratobj.filtered, features = "Percent_ribo", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Post-filtering")
vln.count.patient.after <- 
  VlnPlot(seuratobj.filtered, features = "nCount_RNA", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Post-filtering")
vln.feat.patient.after <- 
  VlnPlot(seuratobj.filtered, features = "nFeature_RNA", pt.size = 0.1,
          split.by = "subtype", cols = subtype_cols) + 
  ggtitle("Post-filtering")

Idents(seuratobj.filtered) <- "subtype"
vln.mt.tumour.after <- 
  customVln(seuratobj.filtered, features = "Percent_mito") + 
  ggtitle("Post-filtering")
vln.rb.tumour.after <- 
  customVln(seuratobj.filtered, features = "Percent_ribo") +
  ggtitle("Post-filtering")
vln.count.tumour.after <- 
  customVln(seuratobj.filtered, features = "nCount_RNA") +
  ggtitle("Post-filtering")
vln.feat.tumour.after <- 
  customVln(seuratobj.filtered, features = "nFeature_RNA") +
  ggtitle("Post-filtering")

Idents(seuratobj.filtered) <- "celltype_major"
vln.mt.cell.after <- 
  customVln(seuratobj.filtered, features = "Percent_mito") + 
  ggtitle("Post-filtering")
vln.rb.cell.after <- 
  customVln(seuratobj.filtered, features = "Percent_ribo") +
  ggtitle("Post-filtering")
vln.count.cell.after <- 
  customVln(seuratobj.filtered, features = "nCount_RNA") +
  ggtitle("Post-filtering")
vln.feat.cell.after <- 
  customVln(seuratobj.filtered, features = "nFeature_RNA") +
  ggtitle("Post-filtering")

# Patchworks of VlnPlots
vln.mt.patient <- (vln.mt.patient.before | vln.mt.patient.after) + 
  plot_annotation(title = "percent.mt by patient ID", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.mt.tumour <- (vln.mt.tumour.before | vln.mt.tumour.after) + 
  plot_annotation(title = "percent.mt by tumour subtype", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.mt.cell <- (vln.mt.cell.before | vln.mt.cell.after) + 
  plot_annotation(title = "percent.mt by major cell type", 
                  theme = title_theme) + 
  plot_layout(guides = "collect")

vln.rb.patient <- (vln.rb.patient.before | vln.rb.patient.after) + 
  plot_annotation(title = "percent.rb by patient ID", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.rb.tumour <- (vln.rb.tumour.before | vln.rb.tumour.after) + 
  plot_annotation(title = "percent.rb by tumour subtype", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.rb.cell <- (vln.rb.cell.before | vln.rb.cell.after) + 
  plot_annotation(title = "percent.rb by major cell type", 
                  theme = title_theme) + 
  plot_layout(guides = "collect")

vln.count.patient <- (vln.count.patient.before | vln.count.patient.after) + 
  plot_annotation(title = "percent.count by patient ID", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.count.tumour <- (vln.count.tumour.before | vln.count.tumour.after) + 
  plot_annotation(title = "percent.count by tumour subtype", 
                  theme = title_theme) + 
  plot_layout(guides = "collect")
vln.count.cell <- (vln.count.cell.before | vln.count.cell.after) + 
  plot_annotation(title = "percent.count by major cell type", 
                  theme = title_theme) + 
  plot_layout(guides = "collect")

vln.feat.patient <- (vln.feat.patient.before | vln.feat.patient.after) + 
  plot_annotation(title = "percent.feat by patient ID", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.feat.tumour <- (vln.feat.tumour.before | vln.feat.tumour.after) + 
  plot_annotation(title = "percent.feat by tumour subtype", 
                  theme = title_theme) + 
  plot_layout(guides = "collect")
vln.feat.cell <- (vln.feat.cell.before | vln.feat.cell.after) + 
  plot_annotation(title = "percent.feat by major cell type", 
                  theme = title_theme) + 
  plot_layout(guides = "collect")

vln.mt <- vln.mt.patient / vln.mt.tumour / vln.mt.cell
vln.rb <- vln.rb.patient / vln.rb.tumour / vln.rb.cell
vln.count <- vln.count.patient / vln.count.tumour / vln.count.cell
vln.feat <- vln.feat.patient / vln.feat.tumour / vln.feat.cell

ggsave(vln.mt, width = 30, height = 25,
       filename = paste0(out.dir, "/plots/violin_percent_mt.pdf"))
ggsave(vln.rb, width = 30, height = 25,
       filename = paste0(out.dir, "/plots/violin_percent_rb.pdf"))
ggsave(vln.count, width = 30, height = 25,
       filename = paste0(out.dir, "/plots/violin_ncount.pdf"))
ggsave(vln.feat, width = 30, height = 25,
       filename = paste0(out.dir, "/plots/violin_nfeature.pdf"))

# Save plots
all.plots <- list(vln.mt.patient.before, vln.rb.patient.before, 
                  vln.count.patient.before, vln.feat.patient.before, 
                  vln.mt.tumour.before, vln.rb.tumour.before, 
                  vln.count.tumour.before, vln.feat.tumour.before, 
                  vln.mt.cell.before, vln.rb.cell.before, 
                  vln.count.cell.before, vln.feat.cell.before,
                  vln.mt.patient.after, vln.rb.patient.after, 
                  vln.count.patient.after, vln.feat.patient.after, 
                  vln.mt.tumour.after, vln.rb.tumour.after, 
                  vln.count.tumour.after, vln.feat.tumour.after, 
                  vln.mt.cell.after, vln.rb.cell.after, 
                  vln.count.cell.after, vln.feat.cell.after)
save(all.plots, file = paste0(out.dir, "/ggplots/seurat_qc.RData"))

# Save filtered cells
write.table(cells.pass, file = snakemake@output[["cells_pass"]],
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

# Save object
saveRDS(seuratobj.filtered, file = snakemake@output[["seurat"]])

# Session info
sessionInfo()