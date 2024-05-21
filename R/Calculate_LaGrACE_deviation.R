#' Scan over a range of penalty value and calculate k-fold LaGrACE deviation for a given dataset.
#'
#' Calculates k-fold LaGrACE deviation over a range of Penalty discount value for optimal penalty
#' discount value selection
#' @param input_data Dataframe containing data for all samples (rows) and variables (columns)
#' @param pd_value Penalty discount value for learning FGES causal network
#' @param maxCat Maximum number of categories for a variable to be treated as discrete
#' @param out_dir Path to directory to which all output files should be written
#' @return A list containing a dataframe of the LaGrACE deviation value (absolute average deviation)
#' and causal graph information (#edges, # connected variables) over a range of Penalty discount value
#' @export


Calculate_LaGrACE_deviation<-function(input_data,pd_value,maxCat=5,out_dir){

  names(input_data) <- make.names(names(input_data),unique = TRUE)

  dir.create(out_dir, showWarnings = FALSE)

  # write out data file with row names
  write.table(input_data, paste0(out_dir,'/FGES_all_samples.w_rownames.txt'), quote = F, sep = '\t', row.names = T, col.names = T)

  # write out data file without row names
  write.table(input_data, paste0(out_dir,'/FGES_all_samples.wout_rownames.txt'), quote = F, sep = '\t', row.names = F, col.names = T)

  PD_Deviation<-c()

  #Perform 10 fold cross validation
  set.seed(123);folds<-caret::createFolds(1:nrow(input_data), k = 10, list = TRUE, returnTrain = T)


  for(i in 1:length(folds)){

    trainData <- input_data[folds[[i]], ]

    # write out reference data file without row names for fges
    dir.create(paste0(out_dir,'/FGES'), showWarnings = FALSE)
    fges_in_file_path <- paste0(paste0(out_dir,'/FGES'),'/train_data_',toString(i),'.txt')
    write.table(trainData, fges_in_file_path, quote = F, sep = '\t', row.names = F, col.names = T)

    lagrace_features_collection<-list()

    for (k in 1:length(pd_value)){
      fges_pd<-pd_value[k]
      # Run FGES
      fges_out_net_file <- learn_fges_network_crossvalidation(fges_in_file_path, fges_pd,maxCat, out_dir=paste0(out_dir,'/FGES'),trainset=i)
      # calculate LaGrACE features from FGES network output
      lagrace_features <- calculate_features(fges_out_net_file, input_data, reference_samples = folds[[i]],maxCat=maxCat)
      lagrace_features$index<-folds[[i]]
      lagrace_features$fges_PD<-fges_pd
      lagrace_features$edges<-nrow(read.table(fges_out_net_file, quote = "", skip = 4, header = F))
      lagrace_features_collection[[k]]<-lagrace_features
      }

    saveRDS(lagrace_features_collection,paste0(out_dir,'/fges_trainset_',toString(i),'_lagrace_features_collection.rds'))

    Mean_Deviation<-matrix(0,nrow = 9,ncol = length(lagrace_features_collection))
    rownames(Mean_Deviation)<-c('# connected feature','abs','abs_connected','abs_unconnected',
                                '# edges','# continuous_connected feature', 'abs_continuous_connected','abs_continuous_unconnected',
                                'FGES_PD')
    for(j in 1:length(lagrace_features_collection)){
      E<-table(lagrace_features_collection[[j]]$variables_selected$num_var!=0)
      Mean_Deviation[1,j]<-E[match("TRUE",names(E))] # number of connected feature
      Mean_Deviation[5,j]<-lagrace_features_collection[[j]]$edges
      Mean_Deviation[9,j]<-lagrace_features_collection[[j]]$fges_PD
      deviation<-lagrace_features_collection[[j]]$features
      index<-setdiff(1:nrow(deviation),lagrace_features_collection[[j]]$index)
      test_deviation<-deviation[index,]
      test_deviation_connected<-test_deviation[,lagrace_features_collection[[j]]$variables_selected$num_var!=0]
      test_deviation_unconnected<-test_deviation[,lagrace_features_collection[[j]]$variables_selected$num_var==0]

      Mean_Deviation[2,j]<-mean(sapply(abs(test_deviation), mean, na.rm = T))
      Mean_Deviation[3,j]<-mean(sapply(abs(test_deviation_connected), mean, na.rm = T))
      Mean_Deviation[4,j]<-mean(sapply(abs(test_deviation_unconnected), mean, na.rm = T))

      connected_continuous_index<-match(lagrace_features_collection[[j]]$continuous_features,names(test_deviation_connected))
      connected_continuous_index<-connected_continuous_index[!is.na(connected_continuous_index)]
      Mean_Deviation[6,j]<-length(connected_continuous_index)

      test_deviation_connected_continuous<-test_deviation_connected[,connected_continuous_index]
      Mean_Deviation[7,j]<-mean(sapply(abs(test_deviation_connected_continuous), mean, na.rm = T))

      unconnected_continuous_index<-match(lagrace_features_collection[[j]]$continuous_features,names(test_deviation_unconnected))
      unconnected_continuous_index<-unconnected_continuous_index[!is.na(unconnected_continuous_index)]
      if (length(unconnected_continuous_index)>0){
        test_deviation_unconnected_continuous<-test_deviation_unconnected[,unconnected_continuous_index]
        Mean_Deviation[8,j]<-mean(sapply(abs(test_deviation_unconnected_continuous), mean, na.rm = T))
      } else {
        test_deviation_unconnected_continuous<-0
        Mean_Deviation[8,j]<-NA
      }
    }
    PD_Deviation<-cbind(PD_Deviation,Mean_Deviation)
  }
  return(PD_Deviation)
}



