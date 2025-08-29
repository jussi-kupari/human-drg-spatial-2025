# === Create a base Seurat object for human nonneurons from Bhuiyan et al. ===
process_bhuiyan_human_nonneurons <- function(seurat_rds) {
  
  # Load data
  seurat <- readRDS(seurat_rds)
  
  DefaultAssay(seurat) <- "RNA"
  seurat@assays$sketch <- NULL
  seurat@assays$integrated <- NULL
  
  # Collect metadata and extract counts
  metadata <- seurat@meta.data
  seurat_data <- seurat@assays[["RNA"]]@counts
  
  # Extract human metadata
  human_metadata <-
    metadata %>% 
    filter(Species == "Human")
  
  # Human cells
  human_cells <- colnames(seurat_data)[colnames(seurat_data) %in% rownames(human_metadata)]
  
  seurat_data <- seurat_data[, human_cells]
  
  
  new_gene_names <- make.unique(rownames(seurat_data))
  rownames(seurat_data) <- new_gene_names
  
  # Create Seurat
  seurat_human <- CreateSeuratObject(seurat_data, meta.data = human_metadata)
  
  
  # Preprocess to PCA
  seurat_human <-
    NormalizeData(seurat_human) %>%
    FindVariableFeatures() %>% 
    ScaleData() %>%
    RunPCA(features = VariableFeatures(seurat_human),
           npcs = 100,
           verbose = FALSE) 
  
  # Run harmony
  seurat_human <- harmony::RunHarmony(seurat_human, "Dataset", max_iter = 50)
  
  # Cluster
  seurat_human <-
    RunUMAP(seurat_human, reduction = "harmony", dims = 1:20) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.5)
  
  seurat_human
}



# === Create a classifier using nonneuronal and immune cells from Su et al.
create_classifier_su_nonneurons <- function(seurat1_rds, seurat2_rds) {
  
  # Load data
  su <- readRDS(seurat1_rds)
  su_immune <- readRDS(seurat2_rds)
  
  # Keep only control cells
  su_nonneurons <- 
    su %>% 
    subset(idents = c("ImmuneCell", "Neuron"), invert = TRUE) %>% 
    SetIdent(value = "treatment") %>% 
    subset(idents = "Control")
  
  su_nonneurons@meta.data <- 
    su_nonneurons@meta.data %>% 
    mutate(celltype = cell_class)
  
  su_immune@meta.data <- 
    su_immune@meta.data %>% 
    mutate(celltype = predicted.id)
  
  su_immune <-
    su_immune %>% 
    SetIdent(value = "treatment") %>% 
    subset(idents = "Control")
  
  # Merge objects
  su_full <- merge(su_nonneurons, su_immune, add.cell.ids = c("nn", "im"))
  
  # Process and run harmony integration
  su_full <-
    su_full %>%
    NormalizeData() %>%
    ScaleData() %>% 
    FindVariableFeatures() %>% 
    RunPCA(npcs = 100)
  
  su_full <- harmony::RunHarmony(su_full, "orig.ident", max_iter = 50)
  
  # Cluster cells
  su_full <-
    RunUMAP(su_full, reduction = "harmony", dims = 1:20) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.5)
  
  # Train classifier
  n_celltypes <- su_full %>% Idents() %>% levels() %>% length()
  cl <- parallel::makePSOCKcluster(n_celltypes)
  doParallel::registerDoParallel(cl)
  
  su_full <-
    su_full %>%
    SetIdent(value = "celltype") %>%
    scPred::getFeatureSpace(pvar = "celltype") %>%
    scPred::trainModel(model = "mda", allowParallel = TRUE)
  
  parallel::stopCluster(cl)
  
  su_full
  
}



# === Predict cell labels for Bhuiyan human nonneurons using Su et al. ===
predict_labels_for_bhuiyan_human_nonneurons <- function(query_seurat, classifier_seurat) {
  
  # Classify
  query_seurat <- scPred::scPredict(query_seurat, classifier_seurat)
  
  # Filter out unassigned cells
  query_filtered <- 
    query_seurat %>% 
    SetIdent(value = "scpred_prediction") %>% 
    subset(idents = "unassigned", invert = TRUE)
  
  # Process to PCA
  query_filtered <-
    query_filtered %>%
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(npcs = 100)
  
  # Harmony integrate
  query_filtered <- harmony::RunHarmony(query_filtered, "Dataset", max_iter = 50)
  
  # Cluster again
  query_filtered <-
    RunUMAP(query_filtered, reduction = "harmony", dims = 1:25) %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = 0.5)
  
  query_filtered
}



# === Humanize gene names for Bhuiyan et al. nonneurons ===
humanize_genes_for_bhuiyan_nonneurons <- function(seurat) {
  
  # Get human orthologs for genes
  genes <-
    babelgene::orthologs(genes = rownames(seurat), species = "mouse", human = FALSE) %>% 
    select(human_symbol, symbol)
  
  # Keep matching orthologs
  features <- 
    tibble(symbol = rownames(seurat)) %>% 
    left_join(genes, join_by(symbol)) %>% 
    distinct(symbol, .keep_all = TRUE) %>% 
    drop_na()
  
  # Filter counts for human genes
  mat <- 
    as.matrix(GetAssayData(seurat, assay = "RNA", layer = "counts")) %>% 
    .[rownames(.) %in% features$symbol, ]
  
  rownames(mat) <- features$human_symbol
  rownames(mat) <- make.unique(rownames(mat))
  
  # Create Seurat
  seurat_humanized <- CreateSeuratObject(counts = mat, meta.data = seurat@meta.data)
  
  # Process Seurat
  seurat_humanized <-
    seurat_humanized %>%
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(npcs = 100) %>%
    harmony::RunHarmony("Dataset", max_iter = 50)
  
  # Cluster Seurat
  seurat_humanized <- 
    RunUMAP(
      seurat_humanized,
      reduction = "harmony",
      dims = 1:25
    ) %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = 0.5)
  
  # Set identities
  Idents(seurat_humanized) <- "scpred_prediction"
  
  # Add humanized as assay to original seurat
  seurat@assays[["RNA_humanized"]] <- seurat_humanized@assays$RNA
  
  # Return original seurat with added assay
  seurat
}
