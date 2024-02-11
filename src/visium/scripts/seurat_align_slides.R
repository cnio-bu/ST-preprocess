log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

## SNAKEMAKE I/O ##
seurat.object <- snakemake@input[["seurat"]]
coord.files <- snakemake@input[["coords_aligned"]]
out.dir <- dirname(snakemake@output[[1]])

## SNAKEMAKE PARAMETERS ##
ref.slide <- snakemake@params[["ref"]]

## FUNCTIONS ##
# Reads aligned coordinates from paste and rounds them to the closest integer
# Also creates a column named after "colname" with the spot names
readCoords <- function(x, colname) {
  df <- read.table(x, sep = "\t", header = TRUE, row.names = 1)  %>%
    round(digits = 0)
  df[[colname]] <- rownames(df)
  return(df)
}

## CODE ##
# Random seed
set.seed(1)

# Read Seurat object
seuratobj <- readRDS(seurat.object)

if (length(coord.files) > 1) {
  # Slides and target slides
  slides <- str_remove(basename(coord.files), pattern = "_.*$")
  target.slides <- slides[which(slides != ref.slide)]
  # List of coordinate files by target slide
  coord.dir <- unique(dirname(coord.files))
  coord.list <- lapply(target.slides, FUN = function(x) {
    list.files(coord.dir, pattern = x, full.names = TRUE)
  })
  names(coord.list) <- target.slides
  # x, y and z axes of reference slide
  xyz.ref <- seuratobj@images[[ref.slide]]@coordinates %>%
    rownames_to_column("spots") %>%
    select(spots, col, row) %>%
    rename(ref = spots,
           x = col,
           y = row) %>%
    mutate(z = str_remove(ref.slide, pattern = "^slide"))
  # Get x, y and z axes of the other slides
  xyz.list <- lapply(seq_along(coord.list), FUN = function(i) {
    filenames <- basename(coord.list[[i]])
    is.ref <- startsWith(filenames, ref.slide)
    target <- readCoords(coord.list[[i]][!is.ref], colname = "target")
    ref <- readCoords(coord.list[[i]][is.ref], colname = "ref")
    coords <- merge(ref, target, by = c("adj_x", "adj_y")) %>% 
      merge(xyz.ref, by = c("ref")) %>%
      select(target, x, y) %>%
      mutate(z = str_remove(names(coord.list)[i], pattern = "^slide"))
    return(coords)
  })
  # Bind all coordinates
  all.coords <- do.call("rbind", xyz.list) %>%
    rename(ref = target) %>%
    rbind(xyz.ref) %>%
    unique() %>%
    arrange(x, y)
  # If there are target spots that match two reference spots, keep the one with 
  # min(x, y)
  all.coords <- all.coords[!duplicated(all.coords$ref), ] %>%
    column_to_rownames("ref")
} else {
  all.coords <- seuratobj@images[[ref.slide]]@coordinates %>%
    select(col, row) %>%
    rename(x = col,
           y = row) %>%
    mutate(z = str_remove(ref.slide, pattern = "^slide"))
}

# Add metadata
seuratobj.aligned <- AddMetaData(seuratobj, metadata = all.coords)

# Save object
saveRDS(seuratobj.aligned, file = snakemake@output[[1]])

# Session info
sessionInfo()