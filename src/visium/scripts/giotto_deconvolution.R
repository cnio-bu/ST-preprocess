log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(Giotto))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
reference <- snakemake@input[["ref"]] 
giotto.object <- snakemake@input[["giotto"]]
out.dir <- dirname(snakemake@output[[1]])

## SNAKEMAKE PARAMETERS ##
n.dim <- snakemake@params[["giotto_ndim"]]

## CODE ##
# Random seed
set.seed(1)

# Patient
patient <- basename(giotto.object) |>
  str_remove(pattern = "_.*$")

# Create outdirs
dir.create(paste0(out.dir, "/plots"), recursive = TRUE, 
                  showWarnings = FALSE)

# Read single-cell reference
subtype <- basename(dirname(reference))
seuratobj <- reference %>%
  str_remove(pattern = "spacexr.*$") %>%
  paste0(., "seurat/qc/subtype/", subtype, "/", subtype, 
         "_single_cell_filtered.rds") %>%
  readRDS()

# Convert to Giotto object
sc.giottoobj <- seuratToGiottoV4(seuratobj, spatial_assay = "RNA")

# Read spatial Giotto object
st.giottoobj <- readRDS(giotto.object)

# Spatial Giotto clustering
st.giottoobj <- calculateHVF(gobject = st.giottoobj)
gene.metadata <- fDataDT(st.giottoobj)
featgenes <- gene.metadata[hvf == "yes"]$feat_ID

st.giottoobj <- runPCA(gobject = st.giottoobj, feats_to_use = featgenes)
elbow <- screePlot(st.giottoobj, ncp = 50)
ggsave(elbow, filename = paste0(out.dir, "/plots/elbow.pdf"))

st.giottoobj <- runUMAP(st.giottoobj, dimensions_to_use = 1:n.dim)
st.giottoobj <- 
  createNearestNetwork(gobject = st.giottoobj, dimensions_to_use = 1:n.dim, k = 15)
st.giottoobj <- 
  doLeidenCluster(gobject = st.giottoobj, resolution = 0.4, n_iterations = 1000)

umap.spatial <- 
  spatDimPlot(gobject = st.giottoobj, cell_color = "leiden_clus", 
              dim_point_size = 2, spat_point_size = 2.5)
ggsave(umap.spatial, 
       filename = paste0(out.dir, "/plots/umap_spatial.pdf"))

# Normalise single-cell counts using the standard method, z-scoring feats over cells
sc.giottoobj <- 
  normalizeGiotto(gobject = sc.giottoobj, scalefactor = 6000, verbose = TRUE)

# Calculate deconvolution signature
markers <- 
  findMarkers_one_vs_all(gobject = sc.giottoobj, method = "scran", 
                         expression_values = "normalized", 
                         cluster_column = "celltype_major")
signature <- unique(markers$feats[which(markers$ranking <= 100)])

# Calculate the mean expression value of signature genes in each cell type
norm.counts <- 
    getExpression(sc.giottoobj, values = "normalized", output = "matrix", 
                  set_defaults = TRUE) |>
    as.matrix()
norm.counts <- 2^(norm.counts[signature, ]) - 1
celltype <- pDataDT(sc.giottoobj)$celltype_major
split.matrix <- split(seq_len(ncol(norm.counts)), celltype)

expr.signature <- 
  lapply(seq_along(split.matrix), function(i) {
    avg <- rowMeans(norm.counts[, split.matrix[[i]], drop = FALSE])
    df <- data.frame(avg)
    colnames(df) <- names(split.matrix)[i]
    return(df)
   }) |>
   bind_cols()

# Run spatialDWLS
st.giottoobj <- 
  runDWLSDeconv(gobject = st.giottoobj, sign_matrix = expr.signature, 
                return_gobject = TRUE)

# Cell type proportions for each spot
cell.prop <- 
  getSpatialEnrichment(st.giottoobj, name = "DWLS", output = "data.table") |>
  as.data.frame() %>%
  column_to_rownames("cell_ID")

# Save cell proportions
write.table(cell.prop, file = snakemake@output[[1]],
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Session info
sessionInfo()