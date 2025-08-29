jie <- readRDS("drg_ra_full.rds")


jie_immune <- readRDS("immune_cells_20221005.rds")



jie_nonneurons <- 
  jie %>% 
  subset(idents = c("ImmuneCell", "Neuron"), invert = TRUE) %>% 
  SetIdent(value = "treatment") %>% 
  subset(idents = "Control")

jie_nonneurons@meta.data <- 
  jie_nonneurons@meta.data %>% 
  mutate(celltype = cell_class)

jie_immune@meta.data <- 
  jie_immune@meta.data %>% 
  mutate(celltype = predicted.id)

jie_immune <-
  jie_immune %>% 
  SetIdent(value = "treatment") %>% 
  subset(idents = "Control")



jie_full <- merge(jie_nonneurons, jie_immune, add.cell.ids = c("nn", "im"))

jie_full <-
  jie_full %>%
  NormalizeData() %>%
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA(npcs = 100)

jie_full <- harmony::RunHarmony(jie_full, "orig.ident", max_iter = 50)

ElbowPlot(jie_full, ndims = 100, reduction = "harmony")

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

jie_full <- cluster_cells(jie_full)

Idents(jie_full) <- "celltype"

DimPlot(jie_full, label = TRUE)

# Train classifier
cl <- parallel::makePSOCKcluster(17)
doParallel::registerDoParallel(cl)

jie_full <-
  jie_full %>%
  SetIdent(value = "celltype") %>%
  scPred::getFeatureSpace(pvar = "celltype") %>%
  scPred::trainModel(model = "mda", allowParallel = TRUE)

parallel::stopCluster(cl)


jie_full %>% scPred::get_scpred()
p <- jie_full %>% scPred::plot_probabilities()


seu <- scPred::scPredict(seu, jie_full)


DimPlot(seu, label = TRUE, repel = TRUE, group.by = "scpred_no_rejection")



seu_filtered <- 
  seu %>% 
  SetIdent(value = "scpred_prediction") %>% 
  subset(idents = "unassigned", invert = TRUE)

seu_filtered <-
  seu_filtered %>%
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 100)

seu_filtered <- harmony::RunHarmony(seu_filtered, "Dataset", max_iter = 50)

ElbowPlot(seu_filtered, ndims = 100, reduction = "harmony")

seu_filtered <- 
  RunUMAP(
    seu_filtered,
    reduction = "harmony",
    dims = 1:25
  ) %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = 0.5)


#DimPlot(seu_filtered, label = TRUE, repel = TRUE, group.by = "scpred_prediction") + 
#DimPlot(seu_filtered, group.by = "Dataset", shuffle = TRUE)

(
  DimPlot(
    seu_filtered,
    label = TRUE,
    repel = TRUE,
    group.by = "scpred_prediction"
  ) +
    NoLegend()
) +
  (DimPlot(seu_filtered, group.by = "Dataset", shuffle = TRUE) + NoLegend()) +
  
  
  (DimPlot(seu_filtered, label = TRUE, repel = TRUE, group.by = "Atlas_annotation") + NoLegend())


saveRDS(seu_filtered, file = "seu_filtered.rds")
