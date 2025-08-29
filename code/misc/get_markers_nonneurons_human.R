
markers <- 
  seu_filtered_human %>% 
  FindAllMarkers(only.pos = TRUE)


top_markers <- 
  markers %>% 
  slice_min(p_val_adj, n = 20, with_ties = FALSE, by = cluster)


DoHeatmap(seu_filtered_human %>% ScaleData(features = top_markers$gene), top_markers$gene)


VlnPlot(seu_filtered_human, top_markers$gene, stack = TRUE, flip = TRUE) + NoLegend()

DotPlot(seu_filtered_human, features = unique(top_markers$gene)) + RotatedAxis()


saveRDS(markers, file = "human_nonneuron_markers.rds")
