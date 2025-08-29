library(tidyverse)
library(Seurat)

data <- readRDS("DRG_nonneurons_release.rds")


DefaultAssay(data) <- "RNA"
data@assays$sketch <- NULL
data@assays$integrated <- NULL

metadata <- data@meta.data
data <- data@assays[["RNA"]]@counts

human_cells <-
  metadata %>% 
  filter(Species == "Human")



keep_these <- colnames(data)[colnames(data) %in% rownames(human_cells)]

data <- data[, keep_these]


new_gene_names <- make.unique(rownames(data))
rownames(data) <- new_gene_names

seu <- CreateSeuratObject(data, meta.data = human_cells)


seu <-
  NormalizeData(seu) %>%
  FindVariableFeatures()

seu <-
  ScaleData(seu) %>%
  RunPCA(features = VariableFeatures(seu),
         npcs = 100,
         verbose = FALSE)  


seu <- harmony::RunHarmony(seu, "Dataset", max_iter = 50)

ElbowPlot(seu, ndims = 100, reduction = "harmony")


cluster_cells <- function(seu) {
  
  seu <-
    RunUMAP(
      seu,
      reduction = "harmony",
      dims = 1:20
    ) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.5)
}

seu <- cluster_cells(seu)

DimPlot(seu, group.by = "Dataset", shuffle = TRUE)
DimPlot(seu, label = TRUE, repel = TRUE, group.by = "Atlas_annotation")




