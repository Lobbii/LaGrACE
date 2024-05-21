#' Calcuate the LaGrACE features for a given dataset
#'
#' Calculates the LaGrACE features for a given dataset based on the provided
#' network learned on a group of reference samples.
#' @param net_file Path to TXT file containing the reference sample network
#' @param input_data Full dataset of with all samples
#' @param reference_samples Logical vector indicating which samples belong to
#' reference group (TRUE) and which do not (FALSE)
#' @param maxCat Maximum number of categories for a variable to be treated as discrete
#' @return A list containing a dataframe of the LaGrACE features, a dataframe
#' indicating which variables were selected for the prediction in calculating each feature
#' , and a dataframe including continuous/discrete feature information
#' @importFrom utils read.table write.table
#' @importFrom stats as.formula lm predict var
#' @export

calculate_features<- function(net_file,input_data,reference_samples,maxCat) {
  data_new<-Convert_Catogorical_to_Binary(input_data,maxCat=maxCat)
  data_binary<-data_new$data_new

  # if there is catogoty:  -1
  if (length(grep("-",colnames(data_binary)))>0){
    colnames(data_binary)<-gsub("-","_",colnames(data_binary))
  }

  var_genes <- names(input_data)
  var_genes_binary<-names(data_binary)

  ref_exp <- data_binary[reference_samples,]

  graph <-read_causality_graph(net_file)


  errorfeatures <- matrix()
  length(errorfeatures) <- dim(data_binary)[1]*length(var_genes_binary)
  dim(errorfeatures) <- c(dim(data_binary)[1],length(var_genes_binary))

  featureselection <- matrix()
  length(featureselection) <- length(var_genes)*2
  dim(featureselection) <- c(length(var_genes),2)
  featureselection<- as.data.frame(featureselection)
  row.names(featureselection) <- var_genes
  names(featureselection) <- c("selected_var","num_var")

  for (i in 1:length(var_genes_binary)){
    Y <- ref_exp[,var_genes_binary[i]]
    Y_full <- data_binary[,var_genes_binary[i]]

    mb <- find_markov_blanket_Pattern(graph =graph ,target =var_genes[data_new$assign[i]])
    mb_index<-data_new$assign %in% match(mb,var_genes)
    mb<-colnames(data_binary)[mb_index]


    if (length(mb)>0){
      W <- ref_exp
      W$Y <- Y
      X <- W[,c(mb,"Y")]

      W <- data_binary
      W$Y_full <- Y_full
      X_full <- W[,c(mb,"Y_full")]
      X_full$Y_full <- NULL

      b <- paste(mb, collapse="+")
      f <- paste("Y ~ ",b,sep="")
      f <- as.formula(f)
      model <- lm(f, data=X)
      predictedY <- predict(model, X_full)

      featureselection[data_new$assign[i],1] <- paste(mb,collapse=",")
      featureselection[data_new$assign[i],2] <- length(mb)

    } else {

      model <- lm(Y~1)

      ## Need dummy variable for X_full here to be a placeholder so we get a vector of the right length for predictedY
      W <- data_binary
      W$Y_full <- Y_full
      X_full <- W[,c(var_genes_binary[i],"Y_full")]
      X_full$Y_full <- NULL
      predictedY <- predict(model, X_full)

      featureselection[data_new$assign[i],1] <- ''
      featureselection[data_new$assign[i],2] <- 0
    }
    errorY <- Y_full - predictedY

    errorfeatures[,i] <- errorY
  }

  errorfeatures <- as.data.frame(errorfeatures)
  names(errorfeatures) <- var_genes_binary
  row.names(errorfeatures) <- row.names(input_data)


  out_list <- list("features" = errorfeatures, "variables_selected" = featureselection,'assign'=data_new$assign,
                   'continuous_features'=data_new$continuous_features) #,'edges'=nrow(graph$edges)

  return(out_list)
}
