#' Cluster samples based on LaGrACE features
#'
#' Uses Seurat graph-based clustering and UMAP/t-SNE embeddings for visualization.
#' Optionally excludes reference samples from clustering.
#'
#' @param features LaGrACE features matrix/data frame (samples x features)
#' @param meta_data_df Data frame containing metadata for all samples (same order as `features`)
#' @param fgs_pd Sparsity penalty value used for FGES (used in output filenames only)
#' @param n_pc Number of principal components used for clustering and embeddings
#' @param total_PC Total number of PCs computed before selecting the first `n_pc`
#' @param cluster_resolution Resolution for Seurat::FindClusters (higher -> more clusters)
#' @param out_dir Output directory for plots and artifacts
#' @param cluster_ref_samples Logical flag: if FALSE, exclude reference samples from clustering
#' @param ref_samples Logical vector indicating reference samples (TRUE = reference)
#' @param seed_value Seed of the random number generator
#' @param Label Optional character vector of labels for plotting (length matches samples)
#' @return List with `clustering` data frame and PCA `loadings` matrix
#' @importFrom grDevices dev.off pdf
#' @export


cluster_sample <- function(features, meta_data_df, fgs_pd, n_pc, total_PC = 50, cluster_resolution, out_dir, cluster_ref_samples = FALSE, ref_samples, seed_value = 1, Label = NULL) {

  if (!is.matrix(features) && !is.data.frame(features)) {
    stop("features must be a matrix or data.frame")
  }
  if (!is.data.frame(meta_data_df)) {
    stop("meta_data_df must be a data.frame")
  }
  if (nrow(meta_data_df) != nrow(features)) {
    stop("meta_data_df must have the same number of rows as features")
  }
  if (!is.logical(ref_samples) || length(ref_samples) != nrow(features)) {
    stop("ref_samples must be a logical vector with length equal to nrow(features)")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (identical(cluster_ref_samples, FALSE)) {
    features <- features[ref_samples == FALSE, , drop = FALSE]
    meta_data_df <- meta_data_df[ref_samples == FALSE, , drop = FALSE]
    if (!is.null(Label)) Label <- Label[ref_samples == FALSE]
  }

  s_obj <- Seurat::CreateSeuratObject(counts = t(features))
  s_obj@meta.data <- cbind(s_obj@meta.data, meta_data_df)

  s_obj@assays$RNA$data <- s_obj@assays$RNA$counts
  s_obj@assays$RNA$scale.data <- as.matrix(s_obj@assays$RNA$counts)

  # Run PCA
  s_obj <- Seurat::RunPCA(object = s_obj, features = rownames(s_obj), npcs = total_PC, seed.use = seed_value, verbose = FALSE)

  pdf(file.path(out_dir, paste0("ssNPAfeatures.fges_pd", toString(fgs_pd), ".pcelbowplot.pdf")))
  print(Seurat::ElbowPlot(object = s_obj))
  dev.off()

  s_obj <- Seurat::FindNeighbors(s_obj, reduction = "pca", dims = 1:n_pc, verbose = FALSE)

  s_obj <- Seurat::FindClusters(s_obj, resolution = cluster_resolution, random.seed = seed_value, verbose = FALSE)

  s_obj <- Seurat::RunTSNE(object = s_obj, dims = 1:n_pc, seed.use = seed_value)

  s_obj <- Seurat::RunUMAP(object = s_obj, dims = 1:n_pc, seed.use = seed_value)

  # pdf(paste0(out_dir, "/ssNPAfeatures.fges_pd", toString(fgs_pd), ".pc_dim", toString(n_pc), ".tsneplot.clusters.pdf"))
  # print(Seurat::DimPlot(s_obj, reduction = "umap"))
  # dev.off()

  # s_obj$Cluster <- Seurat::Idents(object = s_obj)

  s_obj$Cluster <- s_obj$seurat_clusters

  meta_data_name <- names(meta_data_df)[1]

  if (!is.null(Label)){
    s_obj <- Seurat::AddMetaData(s_obj, Label, col.name = "Label")
    my_data <- Seurat::FetchData(s_obj, c("Cluster", meta_data_name, "umap_1", "umap_2", "Label"))
  } else {
    my_data <- Seurat::FetchData(s_obj, c("Cluster", meta_data_name, "umap_1", "umap_2"))
  }

  colnames(my_data)[colnames(my_data)=="umap_1"] <- 'UMAP_1'
  colnames(my_data)[colnames(my_data)=="umap_2"] <- 'UMAP_2'

  output_obj <- list(clustering = my_data, loadings = Seurat::Loadings(s_obj, reduction = "pca")[, 1:n_pc])
  return(output_obj)
}
