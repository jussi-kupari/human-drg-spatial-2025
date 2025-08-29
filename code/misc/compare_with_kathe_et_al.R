library(Seurat)
library(tidyverse)
library(magrittr)
library(scCustomize)

# Load and process Kathe
kathe <-
  ReadMtx(
    "kathe_GSE184370/GSE184370_courtine_lumbar-EES-rehab_UMI.mtx.gz",
    cells = "kathe_GSE184370/GSE184370_barcodes.txt.gz",
    features = "kathe_GSE184370/GSE184370_features.txt.gz", 
    feature.column = 1
  ) %>% 
  CreateSeuratObject() %>% 
  NormalizeData() %>% 
  ScaleData() 

kathe_meta <- 
  read_table("kathe_GSE184370/GSE184370_meta.txt.gz") %>% 
  as.data.frame()

rownames(kathe_meta) <- Cells(kathe)
full_meta <- bind_cols(kathe@meta.data, kathe_meta)
kathe@meta.data <- full_meta
Idents(kathe) <- "global_cell_type"
kathe_neurons <- subset(kathe, idents = "Neurons")
Idents(kathe_neurons) <- "label"

#kathe_uninjured <- subset(kathe_neurons, idents = "Uninjured")

kathe_neurons <-
  kathe_neurons %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# Load Mingdong
md <- readRDS("ref_final_clustering.rds")

# Load protein coding genes
protein_coding <- readRDS("protein_coding.rds")


md <- Add_Mito_Ribo(md, species = "Mouse")
kathe_neurons <- Add_Mito_Ribo(kathe_neurons, species = "Mouse")

md <- Add_Cell_Complexity(md)
kathe_neurons <- Add_Cell_Complexity(kathe_neurons)

# All functions contain
p1 <- 
  QC_Plots_Genes(seurat_object = md, group.by = "treatment") +
  geom_hline(yintercept = median(md@meta.data$nFeature_RNA), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0, 12000) +
  ggtitle(str_glue("median: {round(median(md@meta.data$nFeature_RNA))}"))

p2 <- QC_Plots_UMIs(seurat_object = md, group.by = "treatment") +
  geom_hline(yintercept = mean(md@meta.data$nCount_RNA), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0, 65000) +
  ggtitle(str_glue("median: {round(median(md@meta.data$nCount_RNA))}"))

p3 <- QC_Plots_Mito(seurat_object = md, group.by = "treatment") +
  geom_hline(yintercept = mean(md@meta.data$percent_mito), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0, 20) +
  ggtitle(str_glue("median: {round(median(md@meta.data$percent_mito), digits = 1)}"))

p4 <- QC_Plots_Complexity(seurat_object = md, group.by = "treatment") +
  geom_hline(yintercept = mean(md@meta.data$log10GenesPerUMI), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0.8, 1) +
  ggtitle(str_glue("median: {round(median(md@meta.data$log10GenesPerUMI), digits = 2)}"))

p5 <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 4)

p6 <- QC_Plots_Genes(seurat_object = kathe_neurons, group.by = "global_cell_type") +
  geom_hline(yintercept = median(kathe_neurons@meta.data$nFeature_RNA), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5)  +
  ylim(0, 12000) +
  ggtitle(str_glue("median: {round(median(kathe_neurons@meta.data$nFeature_RNA))}"))

p7 <- QC_Plots_UMIs(seurat_object = kathe_neurons, group.by = "global_cell_type") +
  geom_hline(yintercept = mean(kathe_neurons@meta.data$nCount_RNA), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0, 65000) +
  ggtitle(str_glue("median: {round(median(kathe_neurons@meta.data$nCount_RNA))}"))

p8 <- QC_Plots_Mito(seurat_object = kathe_neurons, group.by = "global_cell_type") +
  geom_hline(yintercept = median(kathe_neurons@meta.data$percent_mito), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0, 20) +
  ggtitle(str_glue("median: {round(median(kathe_neurons@meta.data$percent_mito), digits = 1)}"))

p9 <- QC_Plots_Complexity(seurat_object = kathe_neurons, group.by = "global_cell_type") +
  geom_hline(yintercept = median(kathe_neurons@meta.data$log10GenesPerUMI), 
             linewidth = 3, 
             color = "red", 
             alpha = 0.5) +
  ylim(0.8, 1) +
  ggtitle(str_glue("median: {round(median(kathe_neurons@meta.data$log10GenesPerUMI), digits = 2)}"))

p10 <- patchwork::wrap_plots(p6, p7, p8, p9, ncol = 4)

p5 / p10

(md@meta.data %>% 
    ggplot(aes(nFeature_RNA)) +
    geom_histogram() +
    cowplot::theme_cowplot() +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ggtitle("Current study") +
    
    kathe_neurons@meta.data %>% 
    ggplot(aes(nFeature_RNA)) +
    geom_histogram() +
    cowplot::theme_cowplot() +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ggtitle("Kathe et al.")) /
  
  (md@meta.data %>% 
     ggplot(aes(nCount_RNA)) +
     geom_histogram() +
     cowplot::theme_cowplot() +
     scale_y_continuous(expand = expansion(mult = c(0, .1))) +
     ggtitle("Current study") +
     
     kathe_neurons@meta.data %>% 
     ggplot(aes(nCount_RNA)) +
     geom_histogram() +
     cowplot::theme_cowplot() +
     scale_y_continuous(expand = expansion(mult = c(0, .1))) +
     ggtitle("Kathe et al."))

DimPlot(kathe_neurons, label = TRUE, repel = TRUE, group.by = "cell_type") + NoLegend()

rownames(md)[match(protein_coding, rownames(md))]

protein_coding %in% rownames(md) %>% mean()
protein_coding %in% rownames(kathe_neurons) %>% mean()
rownames(kathe_neurons) %in% rownames(md) %>% mean()
rownames(md) %in% rownames(kathe_neurons) %>% mean()


md_prot <- rownames(md)[match(protein_coding, rownames(md))]
kathe_prot <- rownames(kathe_neurons)[match(protein_coding, rownames(kathe_neurons))]

protein_coding %in% md_prot %>% mean()
protein_coding %in% kathe_prot %>% mean()
kathe_prot %in% md_prot %>% mean()
md_prot %in% kathe_prot %>% mean()


(tibble(kathe = kathe_prot, this_study = md_prot %in% kathe_prot) %>% 
    count(this_study) %>% 
    mutate(percentage = n / sum(n) * 100) %>% 
    ggplot(aes(x = "", y = percentage, fill = this_study)) +
    geom_col(position = "fill") +
    labs(x = "", y = "Prop. genes detected") +
    cowplot::theme_cowplot() +
    scale_y_continuous(labels = scales::percent) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(legend.position = "none")) +
  
  
  (tibble(md = md_prot, kathe = kathe_prot %in% md_prot) %>% 
     count(kathe) %>% 
     mutate(percentage = n / sum(n) * 100) %>% 
     ggplot(aes(x = "", y = percentage, fill = kathe)) +
     geom_col(position = "fill") +
     labs(x = "", y = "Prop. genes detected") +
     cowplot::theme_cowplot() +
     scale_y_continuous(labels = scales::percent) +
     scale_y_continuous(expand = expansion(mult = c(0, .1))) +
     theme(legend.position = "none"))






genes_df <- 
  tibble(protein_coding) %>% 
  mutate(This_study = protein_coding %in% md_prot,
         Kathe_et_al = protein_coding %in% kathe_prot) 

genes_df %>%
  mutate(across(c(This_study, Kathe_et_al), ~ ifelse(., "TRUE", "FALSE"))) %>%
  pivot_longer(cols = c(This_study, Kathe_et_al), names_to = "dataset", values_to = "value") %>%
  group_by(dataset, value) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = dataset, y = percentage, fill = value)) +
  geom_col(position = "fill") +
  labs(x = "Dataset", y = "Prop. mouse prot. coding genes detected") +
  cowplot::theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(legend.position = "none") 


predict_labels <- function(query, reference, dims, k.weight = 50) {
  predictions <- TransferData(
    anchorset = FindTransferAnchors(
      reference = reference, 
      query = query,
      dims = dims
    ),
    k.weight = k.weight,
    refdata = reference@active.ident,
    dims = dims
  )
  
  query %<>%
    AddMetaData(metadata = predictions)
  
  return(query)
}




Idents(kathe_neurons) <- "cell_type"
DefaultAssay(kathe_neurons) <- "RNA"
DefaultAssay(md) <- "RNA"

md <- predict_labels(md, kathe_neurons, dims = 1:20)


(DimPlot(md, label = TRUE, repel = TRUE, group.by = "cluster_id") + 
    NoLegend() + 
    ggtitle("Neuron type (this study)")) +
  (DimPlot(md, label = TRUE, repel = TRUE, group.by = "predicted.id") + 
     NoLegend() +
     ggtitle("Predicted Kathe et al. 2022 neuron type"))


score_umap <-  
  (FeaturePlot(md %>% SetIdent(value = "predicted.id"), 
               "prediction.score.max", 
               order = TRUE, 
               cols = c("red", "green"), 
               label = TRUE,
               repel = TRUE) + ggtitle("Max prediction score"))


(DimPlot(md, label = TRUE, repel = TRUE, group.by = "cluster_id") + 
    NoLegend() + 
    ggtitle("Neuron type (this study)")) +
  (DimPlot(md, label = TRUE, repel = TRUE, group.by = "predicted.id") + 
     NoLegend() +
     ggtitle("Predicted Kathe et al. 2022 neuron type")) + score_umap


prediction_table_mat <- 
  md@meta.data %>% 
  select(cluster_id, starts_with("prediction.score"), -prediction.score.max) %>%
  rename_with(\(x) str_remove(x, "prediction.score.")) %>% 
  pivot_longer(-cluster_id, names_to = "kathe_type", values_to = "prediction_score") %>% 
  count(cluster_id, kathe_type, wt = mean(prediction_score), name = "mean_score") %>% 
  pivot_wider(names_from = kathe_type, values_from = mean_score) %>% 
  column_to_rownames("cluster_id") %>% 
  as.matrix()

clust_md <- hclust(dist(prediction_table_mat)) 
clust_kathe <- hclust(dist(t(prediction_table_mat))) 

heatmap_kathe_to_md <- 
  md@meta.data %>% 
  select(cluster_id, starts_with("prediction.score"), -prediction.score.max) %>%
  rename_with(\(x) str_remove(x, "prediction.score.")) %>% 
  pivot_longer(-cluster_id, names_to = "kathe_type", values_to = "prediction_score") %>% 
  count(cluster_id, kathe_type, wt = mean(prediction_score), name = "mean_score") %>% 
  mutate(
    cluster_id = factor(cluster_id, levels = clust_md$labels[clust_md$order]),
    kathe_type = factor(kathe_type, levels = clust_kathe$labels[clust_kathe$order])
  ) %>% 
  ggplot(aes(cluster_id, kathe_type, fill = mean_score)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::scale_fill_viridis_c() +
  labs(x = "Neuron types (this study)", 
       y = "Kathe et al. neuron types", 
       fill = "Mean score",
       title = "Kathe et al. 2022 to this study")


kathe_neurons <- predict_labels(kathe_neurons, md, dims = 1:20)

(DimPlot(kathe_neurons, label = TRUE, repel = TRUE, group.by = "cell_type") + NoLegend()) +
  (DimPlot(kathe_neurons, label = TRUE, repel = TRUE, group.by = "predicted.id") + NoLegend())

(FeaturePlot(kathe_neurons, 
             "prediction.score.max", 
             order = TRUE, 
             cols = c("red", "green"), 
             label = TRUE,
             repel = TRUE)) +
  (FeaturePlot(kathe_neurons %>% SetIdent(value = "predicted.id"), 
               "prediction.score.max", 
               order = TRUE, 
               cols = c("red", "green"), 
               label = TRUE,
               repel = TRUE))

(FeaturePlot(md %>% SetIdent(value = "predicted.id"), 
             "prediction.score.max", 
             order = TRUE, 
             cols = c("red", "green"), 
             label = TRUE,
             repel = TRUE))


prediction_table_mat_kathe <- 
  kathe_neurons@meta.data %>% 
  select(cell_type, starts_with("prediction.score"), -prediction.score.max) %>%
  rename_with(\(x) str_remove(x, "prediction.score.")) %>% 
  pivot_longer(-cell_type, names_to = "md_type", values_to = "prediction_score") %>% 
  count(cell_type, md_type, wt = mean(prediction_score), name = "mean_score") %>% 
  pivot_wider(names_from = md_type, values_from = mean_score) %>% 
  column_to_rownames("cell_type") %>% 
  as.matrix()

clust2_kathe <- hclust(dist(prediction_table_mat_kathe)) 
clust2_md <- hclust(dist(t(prediction_table_mat_kathe))) 

heatmap_md_to_kathe <- 
  kathe_neurons@meta.data %>% 
  select(cell_type, starts_with("prediction.score"), -prediction.score.max) %>%
  rename_with(\(x) str_remove(x, "prediction.score.")) %>% 
  pivot_longer(-cell_type, names_to = "md_type", values_to = "prediction_score") %>% 
  count(cell_type, md_type, wt = mean(prediction_score), name = "mean_score") %>%  
  mutate(
    cell_type = factor(cell_type, levels = clust2_kathe$labels[clust2_kathe$order]),
    md_type = factor(md_type, levels = clust2_md$labels[clust2_md$order])
  ) %>% 
  ggplot(aes(cell_type, md_type, fill = mean_score)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::scale_fill_viridis_c() +
  labs(x = "Kathe et al. neuron types", 
       y = "Neuron types (this study)", 
       fill = "Mean score",
       title = "This study to Kathe et al. 2022")


heatmap_kathe_to_md + heatmap_md_to_kathe
