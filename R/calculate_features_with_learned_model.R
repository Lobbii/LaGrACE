#' Train LaGrACE model using a given reference dataset and calcuate the LaGrACE features for another dataset
#'
#' This function train LaGrACE predictive model by a given reference dataset 
#' and calculate LaGrACE featur (deviation) for another dataset
#' @param net_file Path to TXT file containing the reference sample network
#' @param input_data_1 Dataset for LaGrACE model training
#' @param reference_samples Logical vector indicating which samples belong to
#' reference group (TRUE) and which do not (FALSE)
#' @param input_data_2 Dataset for LaGrACE features calculation
#' @param maxCat Maximum number of categories for a variable to be treated as discrete
#' @return A list containing a dataframe of the LaGrACE features, a dataframe
#' indicating which variables were selected for the prediction in calculating each feature
#' , and a dataframe including continuous/discrete feature information
#' @importFrom utils read.table write.table
#' @importFrom stats as.formula lm predict var
#' @export

calculate_features_with_learned_model<- function(net_file,input_data_1,reference_samples,input_data_2,maxCat) {
  
  data_train<-Convert_Catogorical_to_Binary(input_data_1,maxCat=maxCat)
  data_binary_train<-data_train$data_new
  
  data_test<-Convert_Catogorical_to_Binary(input_data_2,maxCat=maxCat)
  data_binary_test<-data_test$data_new
  
  if(all(data_train$assign  == data_test$assign)){
    colnames(data_binary_test) <- colnames(data_binary_train)
  }
  
  names(data_binary_test) <- make.names(names(data_binary_test))
  names(data_binary_train) <- make.names(names(data_binary_train))
  table(names(data_binary_train) ==names(data_binary_test)  )
  
  # table( names(data_binary_test) == names(data_binary_train))
  # table( colnames(data_binary_test) == colnames(data_binary_train))
  
  # if there is catogoty:  -1 
  if (length(grep("-",colnames(data_binary_train)))>0){
    colnames(data_binary_train)<-gsub("-","_",colnames(data_binary_train))
  }
  
  var_genes <- names(input_data_1)
  var_genes_binary<-names(data_binary_train)
  
  ref_exp <- data_binary_train[reference_samples,]
  
  graph <-read_causality_graph(net_file)
  

  errorfeatures <- matrix()
  length(errorfeatures) <- dim(data_binary_test)[1]*length(var_genes_binary)
  dim(errorfeatures) <- c(dim(data_binary_test)[1],length(var_genes_binary))
  
  featureselection <- matrix()
  length(featureselection) <- length(var_genes)*2
  dim(featureselection) <- c(length(var_genes),2)
  featureselection<- as.data.frame(featureselection)
  row.names(featureselection) <- var_genes
  names(featureselection) <- c("selected_var","num_var")
  
  for (i in 1:length(var_genes_binary)){
    Y <- ref_exp[,var_genes_binary[i]]
    Y_full <- data_binary_test[,var_genes_binary[i]]
    
    mb <- find_markov_blanket_Pattern(graph =graph ,target =var_genes[data_train$assign[i]])
    mb_index<-data_train$assign %in% match(mb,var_genes)
    mb<-colnames(data_binary_train)[mb_index]
    
    
    if (length(mb)>0){
      W <- ref_exp
      W$Y <- Y
      X <- W[,c(mb,"Y")]
      
      W <- data_binary_test
      W$Y_full <- Y_full
      X_full <- W[,c(mb,"Y_full")]
      X_full$Y_full <- NULL
      
      b <- paste(mb, collapse="+")
      f <- paste("Y ~ ",b,sep="")
      f <- as.formula(f)
      model <- lm(f, data=X)
      predictedY <- predict(model, X_full)
      
      featureselection[data_train$assign[i],1] <- paste(mb,collapse=",")
      featureselection[data_train$assign[i],2] <- length(mb)
      
    } else {
      
      model <- lm(Y~1)
      
      ## Need dummy variable for X_full here to be a placeholder so we get a vector of the right length for predictedY
      W <- data_binary_test
      W$Y_full <- Y_full
      X_full <- W[,c(var_genes_binary[i],"Y_full")]
      X_full$Y_full <- NULL
      predictedY <- predict(model, X_full)
      
      featureselection[data_train$assign[i],1] <- ''
      featureselection[data_train$assign[i],2] <- 0
    }
    errorY <- Y_full - predictedY
    
    errorfeatures[,i] <- errorY
  }
  
  errorfeatures <- as.data.frame(errorfeatures)
  names(errorfeatures) <- var_genes_binary
  row.names(errorfeatures) <- row.names(input_data_2)
  
  
  out_list <- list("features" = errorfeatures, "variables_selected" = featureselection,'assign'=data_train$assign,
                   'continuous_features'=data_train$continuous_features) #,'edges'=nrow(graph$edges)
  
  return(out_list)
}
