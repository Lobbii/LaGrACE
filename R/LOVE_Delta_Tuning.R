#' Tuning Delta value for LOVE algorithm
#'
#' @param Y_features n by p data matrix.
#' @param delta_seq grid of constants for tuning delta(cluster anchor tuning)
#' @return a list including:
#' \itemize{
#' \item delta_seq
#' \item AIC score
#' \item BIC score
#' \item LOVE latent factor
#' }
#' @export
LOVE_Delta_Tuning<-function(Y_features, delta_seq){
  BIC_value<-c()
  AIC_value<-c()
  LOVE_list<-c()
  for (delta in 1:length(delta_seq) ){
    Y_LOVE<-LOVE(Y_features, delta = delta_seq[delta], HT = F, lbd = 1, mu = 1,
                 diagonal = FALSE, stand = TRUE, N  = 50, Union = TRUE,
                 verbose = FALSE)

    LOVE_list<-c(LOVE_list,list(Y_LOVE))

    A_M<-Y_LOVE$A
    Love_features<-MASS::ginv(t(A_M)  %*% A_M) %*% t(A_M)%*%t(as.matrix(Y_features))
    reconstructed<-A_M %*% Love_features

    res<-Y_features - t(reconstructed)
    ssr <- sum(res^2)
    logL <- -(ncol(Y_features)/2)*(log(ssr/ncol(Y_features)))-(ncol(Y_features)/2)*(log(2*pi))-(ncol(Y_features)/2)

    BIC = Y_LOVE$K * log(dim(Y_features)[1]) - 2*logL
    AIC = 2*Y_LOVE$K - 2*logL
    print(c('BIC: ',BIC, ' AIC: ', AIC))
    BIC_value<-c(BIC_value,BIC)
    AIC_value<-c(AIC_value,AIC)
  }

  return(list('Delta' =delta_seq, 'AIC' =AIC_value, 'BIC' =BIC_value, 'LOVE_result' = LOVE_list ))
}
