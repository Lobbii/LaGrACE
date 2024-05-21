#' Filter a dataframe based on the variance of the columns.
#'
#' This function filters a dataframe for the columns with the highest variance.
#' @param df Dataframe containing data for all samples (rows) and variables (columns)
#' @param n_var Number of variables to keep for variance filter
#' @return Dataframe containing the data filtered for highest column variance
#' @export
filter_by_variance <- function(df,n_var) {
  if (n_var < dim(df)[2]){
    v <- apply(df, MARGIN = 2, FUN = var)
    data_sorted <- df[, order(v, decreasing = T)]
    data_subset = data_sorted[, 1:n_var]
    return(data_subset)
  } else {
    return(df)
  }
}
