

genes <-
  babelgene::orthologs(genes = rownames(seu_filtered), species = "mouse", human = FALSE) %>% 
  dplyr::select(human_symbol, symbol) 

features <- tibble(symbol = rownames(seu_filtered))



features <- semi_join(features, genes, join_by(symbol))



features <-
  left_join(features, genes, join_by(symbol)) %>% 
  distinct(symbol, .keep_all = TRUE)


mat <- as.matrix(seu_filtered@assays[["RNA"]]@counts)

mat <- mat[rownames(mat) %in% features$symbol, ]

rownames(mat) <- features$human_symbol
#rownames(mat) <- make.unique(rownames(mat))

Mat <- mat[!is.na(rownames(mat)), ]

seu_filtered_human <- CreateSeuratObject(counts = Mat, meta.data = seu_filtered@meta.data)

seu_filtered_human <-
  seu_filtered_human %>%
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 100)

seu_filtered_human <- harmony::RunHarmony(seu_filtered_human, "Dataset", max_iter = 50)

ElbowPlot(seu_filtered_human, ndims = 100, reduction = "harmony")

seu_filtered_human <- 
  RunUMAP(
    seu_filtered_human,
    reduction = "harmony",
    dims = 1:25
  ) %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = 0.5)


Idents(seu_filtered_human) <- "scpred_prediction"


(
  DimPlot(
    seu_filtered_human,
    label = TRUE,
    repel = TRUE,
    group.by = "scpred_prediction"
  ) +
    NoLegend()
) +
  (DimPlot(seu_filtered_human, group.by = "Dataset", shuffle = TRUE) + NoLegend()) +
  
  
  (DimPlot(seu_filtered_human, label = TRUE, repel = TRUE, group.by = "Atlas_annotation") + NoLegend())

saveRDS(seu_filtered_human, file = "seu_filtered_human.rds")
