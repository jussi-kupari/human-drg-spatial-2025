library(Seurat)
library(tidyverse)
library(magrittr)
library(scPred)

# Load objects

kathe_neurons <- readRDS("kathe_neurons.rds") 
md <- readRDS("ref_final_clustering.rds")

kathe_neurons <- 
  kathe_neurons %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

md <- 
  md %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

md <- getFeatureSpace(md, "cluster_id")
md <- trainModel(md, model = "mda", allowParallel = TRUE)

saveRDS(md, file = "md_with_mda_model.rds")

kathe_neurons <- getFeatureSpace(kathe_neurons, "cell_type")
kathe_neurons <- trainModel(kathe_neurons, model = "mda", allowParallel = TRUE)

saveRDS(kathe_neurons, file = "kathe_neurons_with_mda_model.rds")

stopCluster(cl)