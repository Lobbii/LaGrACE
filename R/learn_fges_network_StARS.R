#' Learn a fast greedy equivalence search (FGES) network with
#' Stability Approach to Regularization Selection (StARS)
#'
#' This function returns the path to a file containing the FGES network learned with StARS
#' on the input dataset.
#' @param fges_in_file_path Path to file containing dataset on which to run FGES
#' @param maxCat Maximum number of categories for a variable to be treated as discrete
#' @param out_dir to directory to which all output files should be written
#' @return Path to file containing the output from FGES
#' @export





learn_fges_network_StARS <- function(fges_in_file_path,maxCat =5, out_dir){
  ssnpa_path <- find.package("LaGrACE")
  s1 <- "java -jar "
  s1_5 <- "/java/tetrad.jar -alg FGES  -d "
  s3 <- "  -stars 10 -highStars 11 -lowStars 1 "
  s4.5 <- " -maxCat "
  s5 <- "  -o "
  s7 <- paste0(out_dir,"/fges_stars.txt")

  fges_java_call <- paste0(s1, ssnpa_path, s1_5, fges_in_file_path, s3,s4.5,toString(maxCat), s5, s7)

  system(fges_java_call)

  fges_net_file <- s7

  return(fges_net_file)
}

