suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(BiocParallel)
  library(progressr)
  library(sf)
})

# helper function to return Seurat object metadata
callmeta <- function (object = NULL) {
  return(object@meta.data)
}

# set directory to read from
dir_use <- "data/"

vizgen.obj <-
  LoadVizgen(
    data.dir = dir_use,
    fov = "merfish.test",
    assay = "Vizgen",
    metadata = c("volume", "fov"),
    # add cell volume info
    type = c("segmentations", "centroids"),
    # type of cell spatial coord matrices
    z = 3L,
    add.zIndex = TRUE,
    # add z slice section to a cell
    update.object = TRUE,
    use.BiocParallel = TRUE,
    workers.MulticoreParam = 14,
    # default 14, for `BiocParallel` processing
    min.area = 5,
    # minimal polygon area to use as a threshold for filtering segmentaion geometries
    add.molecules = TRUE,
    # if to add "molecules" coordinates to FOV of the object
    verbose = TRUE
  )

saveRDS(vizgen.obj, file = "human_drg_vizgen_obj.rds")

vizgen.obj <-
  NormalizeData(vizgen.obj) %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(vizgen.obj)) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:20)

VlnPlot(vizgen.obj, "nCount_Vizgen")


seu <- subset(vizgen.obj, nCount_Vizgen > 400)

seu <-
  NormalizeData(seu) %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(vizgen.obj)) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)



ImageDimPlot(seu, cols = "polychrome", axes = TRUE)
ImageDimPlot(seu, cols = "red", cells = WhichCells(seu, idents = 2))
ImageFeaturePlot(vizgen.obj, features = "SLC17A6")

p1 <- ImageFeaturePlot(seu, features = "SLC17A6")
p2 <- ImageDimPlot(seu, molecules = "SLC17A6", nmols = 10000, alpha = 0.3, mols.cols = "red")
p1 + p2

ImageFeaturePlot(seu, features = "TRPM8", mols.size = 0.5, nmols = 10000)
ImageDimPlot(seu, molecules = "SLC17A6", alpha = 0.3, mols.size = 0.5, mols.cols = "red") +
  ImageDimPlot(seu, molecules = "SOX10", alpha = 0.3, mols.size = 0.5, mols.cols = "green")

ImageDimPlot(seu, 
             alpha = 0.2, 
             mols.size = 0.1, 
             nmols = 10000, 
             molecules = c("SLC17A6"), mols.cols = "red") +

ImageDimPlot(seu, 
             alpha = 0.2, 
             mols.size = 0.1, 
             nmols = 10000, 
             molecules = c("SST"), mols.cols = "blue") +
  
  ImageDimPlot(seu, 
               alpha = 0.2, 
               mols.size = 0.1, 
               nmols = 10000, 
               molecules = c("PVALB"), mols.cols = "green")
  


ImageDimPlot(seu, 
             alpha = 0.1, 
             mols.size = 0.1, 
             nmols = 10000, 
             molecules = c("SLC17A6", "TRPM8", "SST")) +

ImageDimPlot(seu, 
             alpha = 0.1, 
             mols.size = 0.1, 
             nmols = 10000, 
             molecules = c("SLC17A6", "GFRA2", "PVALB"))

saveRDS(seu, file = "human_drg_filtered.rds")
