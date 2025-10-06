#' Calcuate the LaGrACE features for a given dataset
#'
#' Calculates the LaGrACE features for a given dataset based on the provided
#' FGES network learned on a group of reference samples.
#'
#' @param net_file Path to TXT file containing the reference sample network
#' @param input_data Data frame or matrix with samples in rows and variables in columns
#' @param reference_samples Logical vector indicating which samples belong to the
#'   reference group (TRUE) and which do not (FALSE)
#' @param maxCat Maximum number of categories for a variable to be treated as discrete
#' @return A list with:
#'   \itemize{
#'     \item features: data frame of LaGrACE features (residuals)
#'     \item variables_selected: data frame with selected variables and counts per target
#'     \item assign: integer vector that maps binary-expanded columns back to original variables
#'     \item continuous_features: character vector of variables treated as continuous
#'   }
#' @importFrom utils read.table write.table
#' @importFrom stats as.formula lm predict var
#' @export

calculate_features <- function(net_file, input_data, reference_samples, maxCat) {
  # Validate inputs
  if (!is.character(net_file) || length(net_file) != 1L || !file.exists(net_file)) {
    stop("net_file must be an existing path to a TXT network file")
  }
  if (!is.data.frame(input_data) && !is.matrix(input_data)) {
    stop("input_data must be a data.frame or matrix")
  }
  if (!is.logical(reference_samples) || length(reference_samples) != nrow(input_data)) {
    stop("reference_samples must be a logical vector with length equal to nrow(input_data)")
  }
  if (!is.numeric(maxCat) || length(maxCat) != 1L || maxCat < 2) {
    stop("maxCat must be a numeric scalar >= 2")
  }

  # Convert categorical variables to binary indicators
  data_new <- Convert_Catogorical_to_Binary(input_data, maxCat = maxCat)
  data_binary <- data_new$data_new

  # Replace '-' in column names to keep valid model terms
  if (any(grepl("-", colnames(data_binary), fixed = TRUE))) {
    colnames(data_binary) <- gsub("-", "_", colnames(data_binary), fixed = TRUE)
  }

  var_genes <- colnames(input_data)
  var_genes_binary <- colnames(data_binary)

  ref_exp <- data_binary[reference_samples, , drop = FALSE]
  graph <- read_causality_graph(net_file)

  # Preallocate output containers for performance and clarity
  num_samples <- nrow(data_binary)
  num_binary_vars <- length(var_genes_binary)
  errorfeatures <- matrix(0, nrow = num_samples, ncol = num_binary_vars)

  featureselection <- data.frame(
    selected_var = rep("", length(var_genes)),
    num_var = integer(length(var_genes)),
    row.names = var_genes,
    stringsAsFactors = FALSE
  )

  # Fit per-target regressions using Markov blanket predictors
  for (i in seq_along(var_genes_binary)) {
    target_bin_name <- var_genes_binary[i]
    Y <- ref_exp[[target_bin_name]]
    Y_full <- data_binary[[target_bin_name]]

    original_var <- var_genes[data_new$assign[i]]
    mb <- find_markov_blanket_Pattern(graph = graph, target = original_var)
    mb_index <- data_new$assign %in% match(mb, var_genes)
    mb_vars <- colnames(data_binary)[mb_index]

    if (length(mb_vars) > 0) {
      X <- cbind(ref_exp[, mb_vars, drop = FALSE], Y = Y)
      X_full <- data_binary[, mb_vars, drop = FALSE]

      formula_str <- paste("Y ~", paste(mb_vars, collapse = "+"))
      model <- lm(stats::as.formula(formula_str), data = X)
      predictedY <- stats::predict(model, newdata = X_full)

      featureselection[original_var, "selected_var"] <- paste(mb_vars, collapse = ",")
      featureselection[original_var, "num_var"] <- length(mb_vars)
    } else {
      # Intercept-only model; build placeholder design to get correct length
      model <- lm(Y ~ 1)
      X_full <- data_binary[, c(target_bin_name), drop = FALSE]
      predictedY <- stats::predict(model, newdata = X_full)

      featureselection[original_var, "selected_var"] <- ""
      featureselection[original_var, "num_var"] <- 0L
    }

    errorfeatures[, i] <- Y_full - predictedY
  }

  errorfeatures <- as.data.frame(errorfeatures)
  colnames(errorfeatures) <- var_genes_binary
  rownames(errorfeatures) <- rownames(input_data)

  out_list <- list(
    features = errorfeatures,
    variables_selected = featureselection,
    assign = data_new$assign,
    continuous_features = data_new$continuous_features
  )

  return(out_list)
}
