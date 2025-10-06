#' Learn a fast greedy equivalence search (FGES) network
#'
#' This function runs FGES via the bundled Java binary and returns the path to
#' the learned network (TXT format).
#'
#' @param fges_in_file_path Path to TXT file containing the dataset on which to run FGES
#' @param fges_pd Numeric penalty discount value for FGES
#' @param maxCat Integer maximum number of categories for a variable to be treated as discrete
#' @param out_dir Path to directory where FGES outputs will be written
#' @return Path to the FGES network output file
#' @export


learn_fges_network <- function(fges_in_file_path, fges_pd, maxCat = 5, out_dir) {
  # Validate inputs for clearer error messages
  if (!is.character(fges_in_file_path) || length(fges_in_file_path) != 1L) {
    stop("fges_in_file_path must be a single path string")
  }
  if (!file.exists(fges_in_file_path)) {
    stop(sprintf("Input file does not exist: %s", fges_in_file_path))
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

  # Ensure output directory exists
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Locate bundled tetrad jar
  pkg_path <- find.package("LaGrACE")
  jar_primary <- file.path(pkg_path, "java", "tetrad.jar")
  jar_alt <- file.path(pkg_path, "java", "tetrad-censored.jar")
  jar_path <- if (file.exists(jar_primary)) jar_primary else if (file.exists(jar_alt)) jar_alt else NULL
  if (is.null(jar_path)) {
    stop("Unable to locate tetrad jar in package 'LaGrACE/inst/java'.")
  }

  # Build output path
  out_file <- file.path(out_dir, paste0("fges_net_pd_", toString(fges_pd), ".txt"))

  # Execute FGES
  args <- c(
    "-jar", jar_path,
    "-alg", "FGES",
    "-d", fges_in_file_path,
    "-penalty", toString(fges_pd),
    "-maxCat", toString(maxCat),
    "-o", out_file
  )
  res <- tryCatch(
    system2("java", args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) stop(sprintf("Failed to launch FGES: %s", e$message))
  )
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("FGES failed with status %s. Output:\n%s", status, paste(res, collapse = "\n")))
  }

  if (!file.exists(out_file)) {
    stop(sprintf("FGES did not produce expected output file: %s", out_file))
  }

  return(out_file)
}
