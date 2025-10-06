#' Learn FGES network and calculate LaGrACE divergence score
#'
#' This function orchestrates the LaGrACE workflow: it derives a reference-only
#' dataset, runs FGES to learn a causal network on the reference samples, and then
#' computes LaGrACE features for all samples using that network.
#'
#' @param input_data Data frame or matrix with samples in rows and variables in columns
#' @param reference_samples Logical vector indicating which samples belong to the
#'   reference group (TRUE) vs non-reference (FALSE). Must be length `nrow(input_data)`.
#' @param fges_pd Numeric scalar penalty discount for FGES causal network learning
#' @param maxCat Integer scalar: maximum number of categories for a variable to be
#'   treated as discrete when generating binary indicators
#' @param out_dir Path to a directory where all intermediate and output files are written
#' @return A list of LaGrACE features and related metadata (see `calculate_features`)
#' @export


run_LaGrACE <- function(input_data, reference_samples, fges_pd = 1, maxCat = 5, out_dir) {

  # Basic input validation for clearer errors
  if (!is.data.frame(input_data) && !is.matrix(input_data)) {
    stop("input_data must be a data.frame or matrix")
  }
  if (!is.logical(reference_samples) || length(reference_samples) != nrow(input_data)) {
    stop("reference_samples must be a logical vector with length equal to nrow(input_data)")
  }
  if (!is.numeric(fges_pd) || length(fges_pd) != 1L) {
    stop("fges_pd must be a numeric scalar")
  }
  if (!is.numeric(maxCat) || length(maxCat) != 1L || maxCat < 2) {
    stop("maxCat must be a numeric scalar >= 2")
  }
  if (!is.character(out_dir) || length(out_dir) != 1L) {
    stop("out_dir must be a single path string")
  }

  # Ensure unique, syntactically valid column names
  colnames(input_data) <- make.names(colnames(input_data), unique = TRUE)

  # Normalize and create output directories
  out_dir <- sub("/+\$", "", out_dir)  # defensive: strip trailing slashes
  out_dir <- sub("/*$", "", out_dir)    # ensure no trailing slash
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Derive reference-only dataset
  ref_data <- input_data[reference_samples, , drop = FALSE]

  # Write inputs with and without row names for downstream steps
  write.table(
    input_data,
    file = file.path(out_dir, "input_data.all_samples.w_rownames.txt"),
    quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE
  )

  write.table(
    ref_data,
    file = file.path(out_dir, "input_data.reference_samples.w_rownames.txt"),
    quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE
  )

  fges_dir <- file.path(out_dir, "FGES")
  dir.create(fges_dir, showWarnings = FALSE, recursive = TRUE)

  # Write reference data without row names for FGES
  fges_in_file_path <- file.path(fges_dir, "input_data.reference_samples.wout_rownames.txt")
  write.table(
    ref_data,
    file = fges_in_file_path,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  # Run FGES and compute LaGrACE features
  fges_out_net_file <- learn_fges_network(
    fges_in_file_path = fges_in_file_path,
    fges_pd = fges_pd,
    maxCat = maxCat,
    out_dir = fges_dir
  )

  lagrace_features <- calculate_features(
    net_file = fges_out_net_file,
    input_data = input_data,
    reference_samples = reference_samples,
    maxCat = maxCat
  )

  return(lagrace_features)
}
