log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[[1]]

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

# Split object and normalise independently
seuratobjs <- SplitObject(seuratobj, split.by = "orig.ident") |>
  lapply(FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

# Select the features that are repeteadly variable accross samples
features <- SelectIntegrationFeatures(object.list = seuratobjs)

# Find integration anchors and integrate
anchors <- FindIntegrationAnchors(object.list = seuratobjs,
                                  anchor.features = features)
seuratobj.integrated <- IntegrateData(anchorset = anchors)

# Save object
saveRDS(seuratobj.integrated , file = snakemake@output[[1]])

# Session info
sessionInfo()