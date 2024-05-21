setwd('./Simulation_test/Dataset_75_05')
library(permute)
library(MASS)
library(ggfortify)
library(reshape2)
library(ggplot2)
library(LaGrACE)
library(mclust)
library(MetBrewer)
library(infotheo)

cluster_metric<-matrix(nrow = 10,ncol = 6)
colnames(cluster_metric)<-c('LaGrACE_ARI','LaGrACE_AMI','LaGrACE_res','Con_Var_ARI','Con_Var_AMI','Con_Var_res')

set.seed(123)

for (i in seq(0,9)){
  
  # Merge Dataset
  data_0<-read.delim2(paste0('./Data',i,'/data/data0.txt'))
  data_1<-read.delim2(paste0('./Data',i,'/data/data1.txt'))
  data_2<-read.delim2(paste0('./Data',i,'/data/data2.txt'))
  data_3<-read.delim2(paste0('./Data',i,'/data/data3.txt'))
  
  # Generate label
  input_data<-rbind(data_0,data_1,data_2,data_3)
  datainfo<-NULL
  n=1000
  datainfo$labels<-c(rep('Orig',n),rep('case1',n),rep('case2',n),rep('case3',n))
  datainfo<-as.data.frame(datainfo)

  idx<-permute::shuffle(nrow(datainfo))
  input_data<-input_data[idx,]
  datainfo<-datainfo[idx,,drop=FALSE]
  input_data <- as.data.frame(input_data)
  rownames(input_data)<-make.names(rownames(input_data),unique = T)
  
  # Generate PCA plot for continuous variables
  continuous_var<-select_continuous_var(input_data)
  input_data[continuous_var] <- sapply(input_data[continuous_var],as.numeric)
  pca_data <- prcomp((input_data[,continuous_var]), scale. = TRUE)
  p<-autoplot(pca_data, data = datainfo, colour = 'labels')+
    theme_classic()
  cowplot::save_plot(paste0('./Data',i,'/pca_plot_con_Var.png'), p,base_aspect_ratio = 1, base_height = 5)
  
  # pick optimal PD for reference graph
  dir.create(path=paste0('./Data',i,"/Orig"),showWarnings=FALSE)
  dir.create(path=paste0('./Data',i,"/Orig/PD_selection"),showWarnings=FALSE)
  
  reference_samples = datainfo$labels == 'Orig'
  ref_inputdata<-input_data[reference_samples,]
  
  PD_Deviation<-Calculate_LaGrACE_deviation(input_data=ref_inputdata,pd_value=seq(from = 1, to = 10, by=1),maxCat=5,
                                          out_dir=paste0('./Data',i,"/Orig/PD_selection"))
  saveRDS(PD_Deviation,paste0('./Data',i,"/Orig/ref_PD_Deviation.rds"))
  PD_Deviation<-readRDS(paste0('./Data',i,"/Orig/ref_PD_Deviation.rds"))
  
  
  PD_plot<-as.data.frame(t(PD_Deviation[c(9,1,5,2,3),]))
  colnames(PD_plot)[1]<-c('FGES_PD')
  PD_plot$FGES_PD<-as.factor(PD_plot$FGES_PD)
  S<-aggregate(PD_plot[, 2:5], list(PD_plot$FGES_PD), mean)
  opt_PD<-S$Group.1[which.min(S$abs)]
  PD_plot <- melt(PD_plot ,  id.vars = 'FGES_PD', variable.name = 'meanDeviation')
  p<- ggplot(PD_plot, aes(x=FGES_PD,y=value)) + geom_boxplot(aes(fill=meanDeviation)) +
    stat_summary(fun=mean, geom="point", shape=20, size=2, color="gold1", fill="gold1")+
    facet_wrap( ~ meanDeviation, scales="free",ncol=2) + theme_classic()
  cowplot::save_plot(paste0('./Data',i,"/Orig/optimal_PD.png"), p,base_aspect_ratio = 1.5, base_height = 8)
  
  
  ref_index = datainfo$labels == 'Orig'
  LaGrACE_deviation_opt_PD<-run_LaGrACE(input_data, reference_samples = ref_index, fges_pd = opt_PD, out_dir=paste0('./Data',i,"/Orig/Testing"))
  saveRDS(LaGrACE_deviation_opt_PD,paste0('./Data',i,"/Orig/LaGrACE_deviation_opt_PD.rds"))

  
  LaGrACE_feature<-LaGrACE_deviation_opt_PD$features

  ref_index = datainfo$labels == 'Orig'
  seurat_obj<-generate_seuratobj(features=LaGrACE_feature, meta_data_df =datainfo , cluster_ref_samples = FALSE,
                                 ref_samples=ref_index,Label=NULL)
  clustering_tree(seurat_obj, n_pc=20, cluster_resolution =seq(0,1,0.1),out_dir=paste0('./Data',i,"/Orig"),seed_value=1 )
  dir.create(paste0('./Data',i,"/Orig/Clustering"))
  
  # clustering based on LaGrACE with estimated graph  ==========================================
  clu_res<-0.1
  cluster_information<-cluster_sample(features=LaGrACE_feature, meta_data_df = datainfo,
                                      fgs_pd=opt_PD, n_pc=20, cluster_resolution=clu_res,
                                      out_dir=paste0('./Data',i,"/Orig/Clustering"),
                                      ref_samples=ref_index,Label=NULL)
  while (length(unique(cluster_information$clustering$Cluster))<3 || max(table(cluster_information$clustering$Cluster))>1800 ){
    clu_res<-clu_res + 0.1
    cluster_information<-cluster_sample(features=LaGrACE_feature, meta_data_df = datainfo,
                                        fgs_pd=opt_PD, n_pc=20, cluster_resolution=clu_res,
                                        out_dir=paste0('./Data',i,"/Orig/Clustering"),
                                        ref_samples=ref_index,Label=NULL)
  }
  
  LaGrACE_ARI<-aricode::ARI(cluster_information$clustering$Cluste,datainfo$labels[!ref_index])
  print(LaGrACE_ARI)
  
  cluster_metric[(i+1),1]<-LaGrACE_ARI
  cluster_metric[(i+1),2]<-aricode::AMI(cluster_information$clustering$Cluster,datainfo$labels[!ref_index])
  cluster_metric[(i+1),3]<-clu_res
  
  cluster_x_subtype <- table(cluster_information$clustering$Cluster, datainfo$labels[!ref_index])
  cluster_x_subtype
  
  p<-pheatmap::pheatmap(cluster_x_subtype,cluster_rows = F,cluster_cols = F,display_numbers = cluster_x_subtype)
  cowplot::save_plot(paste0('./Data',i,"/Orig/heatmap_LaGrACE.png"), p,base_aspect_ratio = 1, base_height = 4)
  p<-ggplot2::ggplot(cluster_information$clustering,ggplot2::aes(UMAP_1,UMAP_2))+
    ggplot2::geom_point(ggplot2::aes(colour=datainfo$labels[!ref_index]))+
    ggplot2::theme_classic()+
    ggplot2::labs(color="Cluster")+
    scale_fill_manual(values=met.brewer("Stevens", 6))
  cowplot::save_plot(paste0('./Data',i,"/Orig/UMAP_LaGrACE.png"), p,base_aspect_ratio = 1, base_height = 6)
  
  # clustering based on continuous variables  ==========================================
  clu_res<-0.1
  opt_PD<-0
  var_clustering<-cluster_sample(features=input_data[continuous_var], meta_data_df = datainfo,
                                 fgs_pd=opt_PD, n_pc=20, cluster_resolution=clu_res,
                                 out_dir=paste0('./Data',i,"/Orig/Clustering"),
                                 ref_samples=ref_index,Label=NULL)
  while (length(unique(var_clustering$clustering$Cluster))<3 || max(table(var_clustering$clustering$Cluster))>1800 ){
    clu_res<-clu_res + 0.1
    var_clustering<-cluster_sample(features=input_data[continuous_var], meta_data_df = datainfo,
                                   fgs_pd=opt_PD, n_pc=20, cluster_resolution=clu_res,
                                   out_dir=paste0('./Data',i,"/Orig/Clustering"),
                                   ref_samples=ref_index,Label=NULL)
  }
  con_ARI<-aricode::ARI(var_clustering$clustering$Cluste,datainfo$labels[!ref_index])
  print(con_ARI)
  
  cluster_metric[(i+1),4]<-con_ARI
  cluster_metric[(i+1),5]<-aricode::AMI(var_clustering$clustering$Cluste,datainfo$labels[!ref_index])
  cluster_metric[(i+1),6]<-clu_res
  
  cluster_x_subtype <- table(var_clustering$clustering$Cluster, datainfo$labels[!ref_index])
  p<-pheatmap::pheatmap(cluster_x_subtype,cluster_rows = F,cluster_cols = F,display_numbers = cluster_x_subtype)
  cowplot::save_plot(paste0('./Data',i,"/Orig/heatmap_Con_var.png"), p,base_aspect_ratio = 1, base_height = 4)
  p<-ggplot2::ggplot(var_clustering$clustering,ggplot2::aes(UMAP_1,UMAP_2))+
    ggplot2::geom_point(ggplot2::aes(colour=datainfo$labels[!ref_index]))+
    ggplot2::theme_classic()+
    ggplot2::labs(color="Cluster")+
    scale_fill_manual(values=met.brewer("Stevens", 6))
  cowplot::save_plot(paste0('./Data',i,"/Orig/UMAP_con_var.png"), p,base_aspect_ratio = 1, base_height = 6)
  
}


write.csv(cluster_metric,'./cluster_metric.csv',row.names = F)

