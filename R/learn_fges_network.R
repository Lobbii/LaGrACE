#' Learn a fast greedy equivalence search (FGES) network
#'
#' This function returns the path to a file containing the FGES network learned
#' on the input dataset.
#' @param fges_in_file_path Path to file containing dataset on which to run FGES
#' @param fges_pd Penalty discount value for learning FGES causal network
#' @param maxCat Maximum number of categories for a variable to be treated as discrete
#' @param out_dir to directory to which all output files should be written
#' @return Path to file containing the output from FGES
#' @export


learn_fges_network <- function(fges_in_file_path, fges_pd,maxCat =5, out_dir){
  ssnpa_path <- find.package("LaGrACE")
  s1 <- "java -jar "
  s1_5 <- "/java/tetrad.jar -alg FGES  -d "
  s3 <- "  -penalty "
  s4 <- toString(fges_pd)
  s4.5 <- " -maxCat "
  s5 <- "  -o "
  s7 <- paste0("/fges_net_pd_", s4,'.txt')

  fges_java_call <- paste0(s1, ssnpa_path, s1_5, fges_in_file_path, s3, s4,s4.5,toString(maxCat), s5, out_dir, s7)

  system(fges_java_call)

  fges_net_file <- paste0(out_dir,"/fges_net_pd_", s4,'.txt')

  return(fges_net_file)
}
