#' Cluster samples based on LaGrACE feature
#'
#' This function utilize graph-based unsupervised clustering algorithm (Seurat) for sample clustring
#' @param features LaGrACE features matrix/dataframe
#' @param meta_data_df Dataframe containing metadata for all samples
#' @param fgs_pd Sparsity Penalty Value for FGES algorithm
#' @param n_pc number of principal components  for clustering and UMAP embedding
#' @param total_PC total number of principal components for clustering and UMAP embedding
#' @param cluster_resolution Value of the resolution parameter, use a value above (below) 1.0 if
#' you want to obtain a larger (smaller) number of communities.
#' @param out_dir output directory for clustering information and visualization (Dimensional Reduction Plot: UMAP)
#' @param cluster_ref_samples Logical vector (TRUE/FALSE) indicating whether reference samples should be excluded/kept
#' from sample clustering
#' @param ref_samples Logical vector indicating which samples belong to
#' reference group (TRUE) and which do not (FALSE)
#' @param seed_value Seed of the random number generator
#' @param Label a string including label information for Dimensional Reduction Plot(UMAP)
#' @return list including one Dataframe containing clustering information and one Matrix containing PC loading information
#' @importFrom grDevices dev.off pdf
#' @export


cluster_sample <- function(features, meta_data_df, fgs_pd, n_pc,total_PC=50, cluster_resolution, out_dir, cluster_ref_samples=FALSE, ref_samples, seed_value = 1,Label=NULL) {

  if (cluster_ref_samples == F){
    features <- features[ref_samples==F, ]
    meta_data_df <- meta_data_df[ref_samples==F, , drop = F]
    Label <-Label[ref_samples==F]
  }

  s_obj <- Seurat::CreateSeuratObject(counts = t(features))
  s_obj@meta.data<-cbind(s_obj@meta.data,meta_data_df)

  s_obj@assays$RNA$data = s_obj@assays$RNA$counts
  s_obj@assays$RNA$scale.data = as.matrix(s_obj@assays$RNA$counts)

  # Run PCA
  s_obj <- Seurat::RunPCA(object = s_obj, features = rownames(s_obj), npcs = total_PC, seed.use = seed_value, verbose=FALSE)

  pdf(paste0(out_dir, "/ssNPAfeatures.fges_pd", toString(fgs_pd), ".pcelbowplot.pdf"))
  print(Seurat::ElbowPlot(object = s_obj))
  dev.off()

  s_obj <- Seurat::FindNeighbors(s_obj,reduction = "pca", dims = 1:n_pc,verbose = FALSE)

  s_obj <- Seurat::FindClusters(s_obj, resolution = cluster_resolution,random.seed = seed_value,verbose=FALSE)

  s_obj <- Seurat::RunTSNE(object = s_obj, dims = 1:n_pc, seed.use = seed_value)

  s_obj <- Seurat::RunUMAP(object = s_obj, dims = 1:n_pc, seed.use = seed_value)

  # pdf(paste0(out_dir, "/ssNPAfeatures.fges_pd", toString(fgs_pd), ".pc_dim", toString(n_pc), ".tsneplot.clusters.pdf"))
  # print(Seurat::DimPlot(s_obj, reduction = "umap"))
  # dev.off()

  # s_obj$Cluster <- Seurat::Idents(object = s_obj)

  s_obj$Cluster <-s_obj$seurat_clusters

  meta_data_name <- names(meta_data_df)[1]

  if (is.null(Label)==FALSE){
    s_obj<-Seurat::AddMetaData(s_obj,Label,col.name = "Label")
    my_data = Seurat::FetchData(s_obj, c("Cluster", meta_data_name, "umap_1", "umap_2","Label"))
  } else {
    my_data = Seurat::FetchData(s_obj, c("Cluster", meta_data_name, "umap_1", "umap_2"))
  }

  colnames(my_data)[colnames(my_data)=="umap_1"] <- 'UMAP_1'
  colnames(my_data)[colnames(my_data)=="umap_2"] <- 'UMAP_2'

  output_obj <- list("clustering"=my_data, "loadings"=Seurat::Loadings(s_obj, reduction = "pca")[, 1:n_pc])
  return(output_obj)
}
