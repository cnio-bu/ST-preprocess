log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
mtx.dir <- dirname(snakemake@params[["matrix"]])
meta.file <- snakemake@params[["metadata"]]
image.dir <- dirname(snakemake@params[["image"]])
gene.column <- snakemake@params[["gene_column"]]

## CODE ##
# Random seed
set.seed(1)

# Slide and patient
slide <- basename(mtx.dir)
patient <- str_remove(basename(dirname(mtx.dir)), pattern = "_.+$")

# Read 10X matrix
mtx <- Read10X(data.dir = mtx.dir, gene.column = gene.column)
colnames(mtx) <- paste(colnames(mtx), slide, sep = "_")

# Read metadata and pathologist annotation
meta <- read.csv(meta.file, row.names = 1) %>%
  mutate(patient = patient,
         slide = slide)
rownames(meta) <- paste(rownames(meta), slide, sep = "_")

# Filter out artefacts
if ("Classification" %in% colnames(meta)) {
  meta <- meta %>%
  filter(Classification != "Artefact" & !is.na(region))
} else {
  meta <- meta %>%
    mutate(Classification = NA)
}

# Subset matrix and pathologist annotation
intersect.spots <- intersect(colnames(mtx), rownames(meta))
mtx <- mtx[, intersect.spots]
meta <- meta[intersect.spots, , drop = FALSE]

# Create Seurat Spatial object
seuratobj <- CreateSeuratObject(counts = mtx, assay = "Spatial", 
                                meta.data = meta, project = patient)

# Create image
image <- Read10X_Image(image.dir = image.dir,
                       image.name = "tissue_lowres_image.png")
rownames(image@coordinates) <- paste(rownames(image@coordinates), 
                                     slide, sep = "_")
image@key <- paste0(slide, "_")
image <- image[Cells(seuratobj)]
DefaultAssay(image) <- DefaultAssay(seuratobj)

# Add image to Seurat Spatial object
seuratobj@images[[slide]] <- image

# Keep spots in tissue
coordinates <- GetTissueCoordinates(seuratobj, cols = c("imagerow", "imagecol"))
seuratobj <- subset(seuratobj, cells = rownames(coordinates))

# Save object
saveRDS(seuratobj, file = snakemake@output[[1]])

# Session info
sessionInfo()