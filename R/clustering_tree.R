#' Plotting clustering tree
#'
#' This function creates a clustering tree showing the relationship between clusterings at different resolutions.
#' @param seurat_obj Seurat Ojbect
#' @param n_pc number of principal components (total 30 PCs) for clustering
#' @param cluster_resolution scan over cluster resolution (from 0 to 1)
#' @param out_dir output directory for clustering tree plot
#' @param seed_value Seed of the random number generator
#' @import clustree
#' @export

clustering_tree<- function(seurat_obj, n_pc=10, cluster_resolution =seq(0,1,0.1),out_dir='./',seed_value=1 ) {

  # n_pc=10;cluster_resolution =seq(0,1,0.1);out_dir='./trial/';seed_value=1

  seurat_obj <- Seurat::RunPCA(object = seurat_obj, features = rownames(seurat_obj), npcs = 50, seed.use = seed_value, verbose=FALSE)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj,reduction = "pca", dims = 1:n_pc,verbose = FALSE)

  for (resolution in cluster_resolution){
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution,random.seed = seed_value,verbose=FALSE)
  }

  p<-clustree::clustree(seurat_obj, prefix = "RNA_snn_res.")
  cowplot::save_plot(paste0(out_dir,'/clusteree.png'), p,
                     base_aspect_ratio = 0.8, base_height = 10)

  p
}



