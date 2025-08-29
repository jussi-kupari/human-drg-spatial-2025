library(Seurat)
library(tidyverse)
library(scPred)

penn <- readRDS("seurat_final_20230602.rds")

DefaultAssay(penn) <- "RNA"

penn <-
  NormalizeData(penn) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)

# Train classifier
cl <- parallel::makePSOCKcluster(16)
doParallel::registerDoParallel(cl)

penn <-
  penn %>%
  SetIdent(value = "cl.conserv") %>%
  scPred::getFeatureSpace(pvar = "cl.conserv") %>%
  scPred::trainModel(model = "mda", allowParallel = TRUE)

parallel::stopCluster(cl)

saveRDS(penn, file = "penn_with_mda_model.rds")

seu <- scPredict(seu, penn, max.iter.harmony = 50)

Idents(seu) <- "scpred_no_rejection"

ImageDimPlot(seu, cols = "polychrome", axes = TRUE, group.by = "scpred_prediction")
ImageDimPlot(seu, cols = "red", cells = WhichCells(seu, idents = "NP3"))

FindMarkers(seu, ident.1 = "NP2") %>% arrange(desc(avg_log2FC)) %>% View()


ImageDimPlot(seu, 
             alpha = 0.1, 
             mols.size = 0.1, 
             nmols = 10000, 
             molecules = c("SST", "ONECUT1"))
