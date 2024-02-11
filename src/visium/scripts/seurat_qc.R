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
  VlnPlot(x, features = features, pt.size = 0.1, group.by = "slide") + 
    blank_x_theme
}

# Draws custom SpatialFeaturePlot
customSpatialFeaturePlot <- function(x, features) {
  lplots <- SpatialFeaturePlot(x, features = features, pt.size.factor = 8, 
                               combine = FALSE)
  return(lplots)
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
dir.create(paste0(out.dir, "/plots/prefilter"), recursive = TRUE, 
           showWarnings = FALSE)
dir.create(paste0(out.dir, "/plots/postfilter"), recursive = TRUE, 
           showWarnings = FALSE)
dir.create(paste0(out.dir, "/ggplots"), showWarnings = FALSE)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Number of slides and patient
slides <- names(seuratobj@images)
nslides <- length(slides)
patient <- seuratobj@meta.data %>%
  pull(patient) %>%
  unique()

# Percentage mito and ribo
seuratobj[["percent.mt"]] <- 
  PercentageFeatureSet(object = seuratobj, pattern = "^MT-")
seuratobj[["percent.rb"]] <- 
  PercentageFeatureSet(object = seuratobj, pattern = "^RP[SL][[:digit:]]")
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(percent.mt = if_else(is.na(percent.mt), 0, percent.mt),
         percent.rb = if_else(is.na(percent.rb), 0, percent.rb))

# VlnPlots before spot filtering
vln.mt.before <- customVln(seuratobj, features = "percent.mt") + 
  ggtitle("Pre-filtering")
vln.rb.before <- customVln(seuratobj, features = "percent.rb") +
  ggtitle("Pre-filtering")
vln.count.before <- customVln(seuratobj, features = "nCount_Spatial") +
  ggtitle("Pre-filtering")
vln.feat.before <- customVln(seuratobj, features = "nFeature_Spatial") +
  ggtitle("Pre-filtering")

# Number of rows to split the SpatialFeaturePlots into and height of the final
# plot
nrows <- 1
h <- 7*nrows

# SpatialFeaturePlots before spot filtering
sf.mt.before <- customSpatialFeaturePlot(seuratobj, features = "percent.mt")
sf.rb.before <- customSpatialFeaturePlot(seuratobj, features = "percent.rb")
sf.count.before <- customSpatialFeaturePlot(seuratobj, features = "nCount_Spatial")
sf.feat.before <- customSpatialFeaturePlot(seuratobj, features = "nFeature_Spatial")

sf.mt.before <- wrap_plots(sf.mt.before, nrow = nrows)
sf.rb.before <- wrap_plots(sf.rb.before, nrow = nrows)
sf.count.before <- wrap_plots(sf.count.before, nrow = nrows)
sf.feat.before <- wrap_plots(sf.feat.before, nrow = nrows)

ggsave(sf.mt.before, width = 17, height = h,
       filename = paste0(out.dir, "/plots/prefilter/spatial_percent_mt_prefilter.pdf"))
ggsave(sf.rb.before, width = 17, height = h,
       filename = paste0(out.dir, "/plots/prefilter/spatial_percent_rb_prefilter.pdf"))
ggsave(sf.count.before, width = 17, height = h,
       filename = paste0(out.dir, "/plots/prefilter/spatial_ncount_prefilter.pdf"))
ggsave(sf.feat.before, width = 17, height = h,
       filename = paste0(out.dir, "/plots/prefilter/spatial_nfeature_prefilter.pdf"))

# QC variables
qcvars <- c("percent.mt", "percent.rb", "nCount_Spatial", "nFeature_Spatial")

# Set limits
limits.list <- list(mt.limits, rb.limits, count.limits, feat.limits)
limits.list <- lapply(seq_along(limits.list), FUN = function(i) {
  setLimits(limits.list[[i]], feature = qcvars[i], obj = seuratobj)
})
names(limits.list) <- qcvars

# Spot-level filtering
spots <- FetchData(seuratobj, vars = qcvars)
spots.subset <- lapply(names(limits.list), FUN = function(x) {
  if (is.null(limits.list[[x]])) {
    spots.subset <- spots %>% rownames()
  } else {
    spots.subset <- spots %>% 
      filter(get(x) >= limits.list[[x]][1] & get(x) <= limits.list[[x]][2]) %>%
      rownames()
  }
})
spots.pass <- Reduce(intersect, spots.subset)
seuratobj.filtered <- seuratobj[, spots.pass]

# VlnPlots after spot filtering
vln.mt.after <- customVln(seuratobj.filtered, features = "percent.mt") + 
  ggtitle("Post-filtering")
vln.rb.after <- customVln(seuratobj.filtered, features = "percent.rb") +
  ggtitle("Post-filtering")
vln.count.after <- customVln(seuratobj.filtered, features = "nCount_Spatial") +
  ggtitle("Post-filtering")
vln.feat.after <- customVln(seuratobj.filtered, features = "nFeature_Spatial") +
  ggtitle("Post-filtering")

# SpatialFeaturePlots after spot filtering
sf.mt.after <- customSpatialFeaturePlot(seuratobj.filtered, features = "percent.mt")
sf.rb.after <- customSpatialFeaturePlot(seuratobj.filtered, features = "percent.rb")
sf.count.after <- customSpatialFeaturePlot(seuratobj.filtered, features = "nCount_Spatial")
sf.feat.after <- customSpatialFeaturePlot(seuratobj.filtered, features = "nFeature_Spatial")

sf.mt.after <- wrap_plots(sf.mt.after, nrow = nrows)
sf.rb.after <- wrap_plots(sf.rb.after, nrow = nrows)
sf.count.after <- wrap_plots(sf.count.after, nrow = nrows)
sf.feat.after <- wrap_plots(sf.feat.after, nrow = nrows)

ggsave(sf.mt.after, width = 17, height = h,
       filename = paste0(out.dir, "/plots/postfilter/spatial_percent_mt_postfilter.pdf"))
ggsave(sf.rb.after, width = 17, height = h,
       filename = paste0(out.dir, "/plots/postfilter/spatial_percent_rb_postfilter.pdf"))
ggsave(sf.count.after, width = 17, height = h,
       filename = paste0(out.dir, "/plots/postfilter/spatial_ncount_postfilter.pdf"))
ggsave(sf.feat.after, width = 17, height = h,
       filename = paste0(out.dir, "/plots/postfilter/spatial_nfeature_postfilter.pdf"))

# Patchwork of VlnPlos
vln.mt <- (vln.mt.before | vln.mt.after) + 
  plot_annotation(title = "percent.mt", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.rb <- (vln.rb.before | vln.rb.after) +
  plot_annotation(title = "percent.rb", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.count <- (vln.count.before | vln.count.after) +
  plot_annotation(title = "nCount_Spatial", theme = title_theme) + 
  plot_layout(guides = "collect")
vln.feat <- (vln.feat.before | vln.feat.after) +
  plot_annotation(title = "nFeature_Spatial", theme = title_theme) + 
  plot_layout(guides = "collect")

ggsave(vln.mt, width = 14,
       filename = paste0(out.dir, "/plots/violin_percent_mt.pdf"))
ggsave(vln.rb, width = 14,
       filename = paste0(out.dir, "/plots/violin_percent_rb.pdf"))
ggsave(vln.count, width = 14,
       filename = paste0(out.dir, "/plots/violin_ncount.pdf"))
ggsave(vln.feat, width = 14,
       filename = paste0(out.dir, "/plots/violin_nfeature.pdf"))

# Save plots
all.plots <- list(vln.mt.before, vln.rb.before, vln.count.before, 
                  vln.feat.before, sf.mt.before, sf.rb.before, sf.count.before, 
                  sf.feat.before, vln.mt.after, vln.rb.after, vln.count.after, 
                  vln.feat.after, sf.mt.after, sf.rb.after, sf.count.after, 
                  sf.feat.after)
save(all.plots, file = paste0(out.dir, "/ggplots/seurat_qc.RData"))

# Save filtered spots
write.table(spots.pass, file = snakemake@output[["spots_pass"]],
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

# Save object
saveRDS(seuratobj.filtered, file = snakemake@output[["seurat"]])

# Session info
sessionInfo()