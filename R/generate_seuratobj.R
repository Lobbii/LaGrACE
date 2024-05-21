#' Generate Seurat Object
#'
#' This function creates a Seurat object for LaGrACE feature
#' @param features LaGrACE features matrix/dataframe
#' @param meta_data_df Dataframe containing metadata for all samples
#' @param cluster_ref_samples Logical vector (TRUE/FALSE) indicating whether reference samples should be excluded/kept
#' from sample clustering
#' @param ref_samples Logical vector indicating which samples belong to
#' reference group (TRUE) and which do not (FALSE)
#' @param Label a string including label information for Dimensional Reduction Plot(UMAP)
#' @return a Seurat object
#' @export

generate_seuratobj<- function(features, meta_data_df, cluster_ref_samples=FALSE, ref_samples,Label=NULL) {

  if (cluster_ref_samples == F){
    features <- features[ref_samples==F, ]
    meta_data_df <- meta_data_df[ref_samples==F, , drop = F]

    if (is.null(Label)==FALSE) {
      Label <-Label[ref_samples==F]
      }
  }

  s_obj <- Seurat::CreateSeuratObject(counts = t(features))
  s_obj@meta.data<-cbind(s_obj@meta.data,meta_data_df)


  if(is.null(Label)==FALSE){
    s_obj@meta.data<-cbind(s_obj@meta.data,Label)
  }


  s_obj@assays$RNA$data = s_obj@assays$RNA$counts
  s_obj@assays$RNA$scale.data = as.matrix(s_obj@assays$RNA$counts)

  return(s_obj)
}
