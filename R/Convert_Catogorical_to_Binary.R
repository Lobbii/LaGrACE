#' Convert discrete variables into binary indicator variables
#'
#' Transforms discrete variables into one-hot encoded indicators, preserving
#' continuous variables as-is. Returns the transformed data, a vector of
#' continuous variable names, and a mapping from transformed columns back to
#' the corresponding original variable index.
#'
#' @param data Data frame or matrix with samples in rows and features in columns
#' @param maxCat Integer maximum number of categories for a variable to be
#'   treated as discrete
#' @return A list with:
#'   \itemize{
#'     \item data_new: data frame with continuous variables and one-hot encodings
#'     \item continuous_features: character vector of continuous variable names
#'     \item assign: integer vector mapping each column in `data_new` to the original variable index
#'   }
#' @export

Convert_Catogorical_to_Binary <- function(data, maxCat = 5) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix")
  }
  if (!is.numeric(maxCat) || length(maxCat) != 1L || maxCat < 2) {
    stop("maxCat must be a numeric scalar >= 2")
  }

  data_frame <- as.data.frame(data, stringsAsFactors = FALSE)
  transformed <- list()
  continuous_features <- character(0)
  assigns <- integer(0)

  for (i in seq_len(ncol(data_frame))) {
    feature <- data_frame[[i]]
    feature_name <- colnames(data_frame)[i]
    unique_vals <- unique(feature)
    is_discrete <- (length(unique_vals) / length(feature) < 0.05) && (length(unique_vals) <= maxCat)

    if (is_discrete) {
      if (length(unique_vals) > 2) {
        dummies <- as.data.frame(varhandle::to.dummy(feature, prefix = feature_name), stringsAsFactors = FALSE)
        for (col_name in colnames(dummies)) {
          transformed[[col_name]] <- dummies[[col_name]]
        }
        assigns <- c(assigns, rep.int(i, ncol(dummies)))
      } else {
        # Binary category -> convert to 0/1
        binary_feature <- as.numeric(as.factor(feature)) - 1
        transformed[[feature_name]] <- binary_feature
        assigns <- c(assigns, i)
      }
    } else {
      transformed[[feature_name]] <- feature
      continuous_features <- c(continuous_features, feature_name)
      assigns <- c(assigns, i)
    }
  }

  data_new <- as.data.frame(transformed, stringsAsFactors = FALSE)
  return(list(data_new = data_new, continuous_features = continuous_features, assign = assigns))
}
