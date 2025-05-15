# LaGrACE
LaGrACE: Estimating gene program dysregulation using latent gene regulatory network



## About LaGrACE
LaGrACE is a novel approach designed to estimate regulatory network-based pseudo control outcome to characterize gene program dysregulation for samples within treatment (or disease) group. This method enables grouping of samples with similar patterns of gene program dysregulation, thereby facilitating discovery of underlying molecular mechanisms induced by treatment or disease.


## Reference
If you use `LaGrACE` in your work, please cite

> **LaGrACE: Estimating gene program dysregulation with latent regulatory networks**
>
> Minxue Jia, Haiyi Mao, Mengli Zhou, Yu-Chih Chen, Panayiotis(Takis) Benos
>
> _Molecular Systems Biology_ . doi: [10.1038/s44320-025-00115-3](https://doi.org/10.1038/s44320-025-00115-3).



## Installation
R >= 4.3.2, Java >= 1.7.0

```r
devtools::install_github("Lobbii/LaGrACE",build_vignettes = FALSE)
```

## Usage

LaGrACE comprises four main steps, which are summarized below and elaborated in the following sections.

(1)	Infer biologically meaningful latent factors (gene programs) from gene expression profiles of reference samples and then calculate factor matrix for treatment group using loading matrix inferred based on reference samples.

(2)	Learn a Bayesian network, based on inferred latent factors and clinical variables of interest, as the reference graph.  This graph encapsulates interactions among latent factors from gene expression profiles and clinical variables.

(3)	Estimate pseudo control outcome to approximate pre-treatment state on a gene program or a clinical measurement by utilizing its Markov blanket derived from the reference graph.

(4)	Calculate difference between pseudo-control outcome and observed value for treatment samples.


**Step 1: Gene Program Discovery**

we have a matrix Y (samples X fetures) and Sample Lable (Referece or Disease)
```r
library(LaGrACE)
library(clustree)
library(MASS)
library(reshape2)

#determine reference sample index
reference_index<- Label %in% 'Reference'  
Y_ref<-Y[reference_index,]

# perform z score normalization based on reference group
Y_mean<-apply(Y_ref, 2, function(x) mean(x))
Y_sd<-apply(Y_ref, 2, function(x) sd(x))
Y_zscore<-sweep(Y, MARGIN = 2, STATS = Y_mean, `-`)
Y_zscore<-sweep(Y_zscore, MARGIN = 2, STATS = Y_sd, `/`)

# scan over possible cluster anchor tuning hyperparameter, delta 
delta_seq<-seq(0.1,3,by=0.1) 
Y_LOVE_result<-LOVE_Delta_Tuning(Y_zscore[reference_index,], delta_seq) 

# pick up optimum delta value based on AIC score
Y_LOVE_opt<-Y_LOVE_result$LOVE_result[[which.min(Y_LOVE_result$AIC)]] 
```

**Step 2 - 4**
```r
# infer latent factors for all samples
A_M<-Y_LOVE_opt$A
Love_features<-ginv(t(A_M)  %*% A_M) %*% t(A_M)%*%t(as.matrix(Y_zscore)) 


rownames(Love_features)<-paste0('LOVE_',1:nrow(Love_features))
input_feature<-as.data.frame(t(Love_features))
names(input_feature)<-make.names(colnames(input_feature),unique = TRUE)
rownames(input_feature)<-make.names(rownames(input_feature),unique = TRUE)


# graph sparsity parameter selection 
dir.create(path=paste0("./PD_selection"),showWarnings=FALSE)
PD_Deviation<-Calculate_LaGrACE_deviation(input_data=input_feature[reference_index,],
                                        pd_value=seq(from = 1, to = 10, by=0.5), 
                                        maxCat=5,out_dir=paste0("./PD_selection"))

PD_plot<-as.data.frame(t(PD_Deviation[c(9,1,5,2,3),]))
colnames(PD_plot)[1]<-c('FGES_PD')
PD_plot$FGES_PD<-as.factor(PD_plot$FGES_PD)
S<-aggregate(PD_plot[, 2:5], list(PD_plot$FGES_PD), mean)
opt_PD<-S$Group.1[which.min(S$abs)] # PD value with smallest prediction error

# Infer network deviation for all samples
LaGrACE_deviation_opt_PD<-run_LaGrACE(input_feature, reference_samples = reference_index, fges_pd = opt_PD, out_dir="./Testing") 
```




