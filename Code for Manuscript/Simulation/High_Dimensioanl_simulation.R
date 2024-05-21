library(ExtDist)
library(rCausalMGM)
library(dplyr)

setwd('./Simulation_test/Dataset_75_05_HD/')

select_continuous_var<-function(dataset){
  maxCat = 4
  continuous_features<-c()
  for (i in 1:ncol(dataset)){
    feature<-dataset[[i]]
    if (length(unique(feature))/length(feature) > 0.05 & length(unique(feature)) > maxCat){
      continuous_features<-c(continuous_features,i)
    } 
  }
  return(continuous_features)
}


simulate_data<-function(data_subset,contunious_var){
  data_subset_con<-data_subset[,contunious_var]
  data_subset_dis<-data_subset[,-contunious_var]
  data_subset_dis <- sapply(data_subset_dis,as.numeric)
  
  loadings<-rnorm(2500*50)
  loadings[abs(loadings)<quantile(abs(loadings),probs = 0.95)]=0
  hist(loadings)
  loadings<-matrix(loadings,nrow=50)
  
  data_subset_con<-as.matrix(data_subset_con)
  simulated_data<-data_subset_con %*% loadings + rnorm(4000*2500)
  
  simulation_list<-list("continuous_data" =data_subset_con, 
                        "data_subset_dis" = data_subset_dis,
                        "loading" = loadings,
                        "simulated_data" = simulated_data,
                        "contunious_var" =colnames(data_subset)[contunious_var] )
  
  return(simulation_list)
}

for (i in seq(0,9)){
  # Merge Dataset
  data_0<-read.delim2(paste0('./Data',i,'/data/data0.txt'))
  data_1<-read.delim2(paste0('./Data',i,'/data/data1.txt'))
  data_2<-read.delim2(paste0('./Data',i,'/data/data2.txt'))
  data_3<-read.delim2(paste0('./Data',i,'/data/data3.txt'))
  
  n=1000 # sample number
  input_data<-rbind(data_0,data_1,data_2,data_3)
  datainfo<-NULL
  datainfo$labels<-c(rep('Orig',n),rep('case1',n),rep('case2',n),rep('case3',n))
  datainfo<-as.data.frame(datainfo)
  
  input_data <- as.data.frame(input_data)
  rownames(input_data)<-make.names(rownames(input_data),unique = T)
  
  contunious_var<-select_continuous_var(input_data)
  input_data[contunious_var] <- sapply(input_data[contunious_var],as.numeric)
  
  dir.create(path=paste0('./Data',i,"/HD"),showWarnings=FALSE)
  saveRDS(datainfo,paste0('./Data',i,"/HD",'/datainfo.rds'))
  saveRDS(input_data,paste0('./Data',i,"/HD",'/input_data.rds'))
  
  simulated_dataset<-simulate_data(input_data,contunious_var)
  saveRDS(simulated_dataset,paste0('./Data',i,"/HD",'/simulated_dataset.rds'))

}
