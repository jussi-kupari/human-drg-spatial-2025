library(targets)

source(here::here("code", "functions.R"))

tar_option_set(packages = c("Seurat", "tidyverse"))

list(
  # Data files
  tar_target(
    bhuiyan_nonneurons_drg,
    here::here("data", "DRG_nonneurons_release.rds"),
    format = "file"
  ),
  tar_target(drg_ra_full, here::here("data", "drg_ra_full.rds"), format = "file"),
  tar_target(
    drg_ra_immune,
    here::here("data", "immune_cells_20221005.rds"),
    format = "file"
  ),
  tar_target(
    penn_drg_updated,
    here::here("data", "seurat_final_20230602.rds"),
    format = "file"
  ),
  
  # Wrangle Bhuiyan et al. human nonneurons to Seurat
  tar_target(
    bhuiyan_nonneurons_seurat,
    process_bhuiyan_human_nonneurons(bhuiyan_nonneurons_drg)
  ),
  
  # Create classifier from Su et al.
  tar_target(
    classifier_nonneurons_su,
    create_classifier_su_nonneurons(drg_ra_full, drg_ra_immune)
  ),
  
  # Classify cells from Bhuiyan et al.
  tar_target(
    bhuiyan_nonneurons_classified,
    predict_labels_for_bhuiyan_human_nonneurons(bhuiyan_nonneurons_seurat, 
                                                classifier_nonneurons_su) %>%
      # Add RNA assay with humanized genes
      humanize_genes_for_bhuiyan_nonneurons()
  )
  
)