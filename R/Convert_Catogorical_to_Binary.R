#' Convert discrete variable into binary indicator variables in a dataframe
#'
#' This function transform all discrete variables into several binary indicator variables,
#' one for each category of the original discrete variable
#' @param data a matrix or dataframe, each row represents a sample, each column represents a feature
#' @param maxCat Full dataset of with all samples
#' @return A list containing a dataframe of the transformed dataset, a dataframe
#' indicating which variables are continuous, and a dataframe including
#' original variable names of transformed data
#' , and a dataframe including continuous/discrete feature information
#' @export

Convert_Catogorical_to_Binary<-function(data,maxCat=5){
  data_new<-matrix(nrow = nrow(data),ncol = 0)
  data_new<-as.data.frame(data_new,stringsAsFactors = F)
  num_discrete<-0
  continuous_features<-c()
  assigns<-c()
  for (i in 1:ncol(data)){
    feature<-data[[i]]
    if (length(unique(feature))/length(feature) <0.05 & length(unique(feature)) <= maxCat){

      if( length(unique(feature))>2  ){
        dummies<-varhandle::to.dummy(feature,prefix = colnames(data)[i])
        data_new<-cbind(data_new,dummies,stringsAsFactors = FALSE)
        assigns<-c(assigns,rep(i,ncol(dummies)))
        num_discrete<-num_discrete +1
      } else {
        # if feature is discrete variable with 2 categories (character)
        binary_feature<-as.numeric(as.factor(data[[i]])) -1
        data_new<-cbind(data_new,binary_feature,stringsAsFactors=FALSE)
        colnames(data_new)[ncol(data_new)]<-colnames(data)[i]
        assigns<-c(assigns,rep(i,1))
        num_discrete<-num_discrete +1
      }

    } else{
      data_new<-cbind(data_new,data[[i]],stringsAsFactors=FALSE)
      colnames(data_new)[ncol(data_new)]<-colnames(data)[i]
      continuous_features<-c(continuous_features,colnames(data)[i])
      assigns<-c(assigns,rep(i,1))
    }
  }
  return(list('data_new' = as.data.frame(data_new), 'continuous_features' =continuous_features,'assign'= assigns))
}
