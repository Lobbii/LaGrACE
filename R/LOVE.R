########################################################################################
#############              LOVE --- Overlapping Clustering               ###############
########################################################################################

#' Run Love Algorithm
#'
#' Non-overlapping/Overlappig variable clustering for given data X where each row denotes a data point and columns for variables.
#' Reference: Bing X, Bunea F, Royer M, Das J. Latent Model-Based Clustering for Biological Discovery. iScience. 2019 Apr 26;14:125-135. doi: 10.1016/j.isci.2019.03.018. Epub 2019 Mar 21. PMID: 30954780; PMCID: PMC6449745.
#' @param X n by p data matrix.
#' @param delta grid of leading constants for tuning delta(cluster anchor tuning). This grid needs to be multiplied by sqrt(log(p)/n)
#' @param HT If TRUE, then hard-thresholding is used for estimating A; else soft-thresholding is used
#' @param lbd membership matrix tuning
#' @param mu the grid of leading constants for tuning the thresholding level diagonal: TRUE if the cov(Z) is diagonal; else FALSE. If TRUE, OLS is used for solving A
#' @param diagonal  TRUE if the cov(Z) is diagonal; else FALSE. If TRUE, #     OLS is used for solving A
#' @param stand True if centering is needed
#' @param N the number of cross validation used for delta
#' @param Union True if using "or" rule for estimating pure nodes; False if using "and" rule
#' @param verbose if the process needs to be printed out
#' @return a list including:
#' \itemize{
#' \item K the number of clusters.
#' \item A the estimated p by K matrix.
#' \item group: the estimated group structure(indices for each cluster).
#' \item pure the vector of estimated pure variables.
#' }
#' @export
LOVE <- function(X, delta = seq(1.5, 2.5 ,0.1), HT = F, lbd = 1, mu = 1,
                 diagonal = FALSE, stand = TRUE, N = 50, Union = TRUE,
                 verbose = TRUE) {

  if (stand)
    X <- scale(X, TRUE, FALSE)
  n <- nrow(X)
  p <- ncol(X)
  deltaGrids <- delta * sqrt(log(p) / n)
  if (length(deltaGrids) > 1) {
    if (verbose)
      cat("Selecting optimal delta by using data splitting...\n")
    optDelta <- stats::median(replicate(N, CV_Delta(X, deltaGrids, Union, diagonal)))
    print(paste0("use delta: ",optDelta/sqrt(log(p) / n)) )
  }else
    optDelta <- deltaGrids
  if (verbose) {
    cat("Finishing selecting optimal delta =", optDelta, "with leading
         constant", min(which(deltaGrids >= optDelta)),"...\n")
    cat("Estimating pure rows...\n")
  }
  Sigma <- t(X) %*% X / nrow(X)
  resultAI <- EstAI(Sigma, optDelta, Union)
  AI <- resultAI$AI
  if (length(resultAI$pureVec) == nrow(Sigma)) {
    group <- resultAI$pureSignInd
  } else {
    if (verbose)
      cat("Estimating C and Sigma_IJ...\n")
    C_hat <- EstC(Sigma, AI, diagonal)
    while (is.null(C_hat)) {
      optDelta <- optDelta * 1.01
      resultAI <- EstAI(Sigma, optDelta, Union)
      AI <- resultAI$AI
      C_hat <- EstC(Sigma, AI, diagonal)
    }
    Y <- EstY(Sigma, AI, resultAI$pureVec)
    if (diagonal) {
      if (verbose)
        cat("Estimating non-pure rows by using Ordinary Least Square...\n")
      C_inv <- diag(diag(C_hat)^{-1})
      AJ <- threshA(t(C_inv %*% Y), mu * optDelta)
    } else {
      if (verbose)
        cat("Selecting the tuning parameter for estimating Omega...\n")
      optLbd <- lbd * optDelta
      print(paste("lbd",lbd,"optDelta",optDelta))
      if (length(lbd) > 1)
        optLbd <- stats::median(replicate(N, CV_lbd(X, lbd * optDelta, AI,
                                             resultAI$pureVec, diagonal)))
      if (verbose) {
        cat("Selecting the optimal lambda =", optLbd, "...\n")
        cat("Estimating Omega...\n")
      }
      Omega <- EstOmega(optLbd, C_hat)
      if (verbose && HT) {
        cat("Estimating non-pure rows by Hard Thresholding...\n")
      } else if (verbose) {
        cat("Estimating non-pure rows by Soft Thresholding...\n")
      }
      if (HT) {
        AJ = threshA(t(Omega %*% Y), mu * optDelta * norm(Omega, "I"))
      } else {
        AJ = EstAJInv(Omega, Y, mu * optDelta * norm(Omega, "I"))
      }
    }
    AI[-resultAI$pureVec, ] <- AJ
    group <-recoverGroup(AI)
  }
  return(list(K = ncol(AI), A = AI, group = group, pure = resultAI$pureVec,
              C = C_hat))
}


#################################################################################
#####                                                                       #####
##### This is the code for Step 2: regression to estimate the non-pure rows #####
#####                                                                       #####
#################################################################################

#' Estimate C. If diagonal=True, estimate only diagonal elements.
#'
#' @param Sigma p by p covariance matrix.
#' @param AI p by K matrix.
#' @param diagonal TRUE or FALSE.
#' @return K by K estimated C_hat
#' @export
EstC <- function(Sigma, AI, diagonal) {
  K <- ncol(AI)
  C <- diag(0, K, K)
  for (i in 1:K) {
    groupi <- which(AI[ ,i] != 0)
    if (length(groupi) == 1) {
      cat("Singleton pure node is not allowed...Need to reselect delta...\n")
      return(NULL)
    }
    sigmai <- as.matrix(abs(Sigma[groupi,groupi]))
    tmpEntry <- sum(sigmai) - sum(diag(sigmai))
    if (length(groupi) == 1) {
      C[i,i] <- tmpEntry
    } else {
      C[i,i] <- tmpEntry / (length(groupi) * (length(groupi) - 1))
    }
    if (!diagonal && i<K) {
      for (j in (i+1):K) {
        groupj <- which(AI[ ,j]!=0)
        # adjust the sign for each row
        sigmaij <- AI[groupi,i] * as.matrix(Sigma[groupi,groupj])
        sigmaij <- t(AI[groupj,j] * t(sigmaij))
        tmpEntry <- sum(sigmaij)
        C[i,j] <- C[j,i] <- tmpEntry / (length(groupi) * length(groupj))
      }
    }
  }
  return(C)
}



#' Estimate the Sigma_{TJ}_hat
#'
#' @param Sigma p by p matrix.
#' @param AI p by K matrix.
#' @param pureVec the vector of pure node indices.
#' @return Sigma_{TJ}_hat K by |J| matrix.
#' @export
EstY <- function(Sigma, AI, pureVec) {
  SigmaS <- AdjustSign(Sigma, AI)
  SigmaJ <- matrix(SigmaS[ ,-pureVec], nrow <- nrow(Sigma))
  SigmaTJ <- matrix(0, ncol(AI), nrow(AI) - length(pureVec))
  for (i in 1:ncol(AI)) {
    groupi <- which(AI[ ,i] != 0)
    SigmaiJ <- as.matrix(SigmaJ[groupi, ])
    SigmaTJ[i, ] <- apply(SigmaiJ, 2, mean) # Average columns along the rows.
  }
  return(SigmaTJ)
}



#' Sign operation on matrix {@code Sigma} according to the sign of {@code AI}
#'
#' @param Sigma p by p matrix;
#' @param AI p by K matrix
#' @return p by p matrix
#' @export
AdjustSign <- function(Sigma, AI) {
  SigmaS <- matrix(0, nrow(AI), nrow(AI))
  for (i in 1:nrow(AI)) {
    index <- which(AI[i, ] != 0)
    if (length(index) != 0)
      SigmaS[i, ] = sign(AI[i,index]) * Sigma[i, ]
  }
  return(SigmaS)
}

#' This function estimates the |J| by K matrix A_J by using soft thresholding.
#'
#' @param Omega K by K estimated C^{-1}.
#' @param Y K by |J| reponse matrix.
#' @param lbd the chosen constant of the RHS constraint for soft-thresholding.
#' @return |J| by K matrix A_J.
#' @export
EstAJInv <- function(Omega,Y,lbd) {
  AJ <- matrix(0, ncol(Y), nrow(Y))
  for (i in 1:ncol(Y)) {
    Atilde <- Omega %*% as.matrix(Y[ ,i])
    AJ[i, ] <- LP(Atilde, lbd)
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
  }
  return(AJ)
}




################################################################################
######                                                                   #######
######   These are the code which convert the original problem into LP   #######
######                                                                   #######
################################################################################
#' This function solve LP program
#'
#' @param col_ind K by K estimated C^{-1}.
#' @param C K by |J| reponse matrix.
#' @param lbd the chosen constant of the RHS constraint for soft-thresholding.
#' @return |J| by K matrix A_J.
#' @export
Solve_row <- function(col_ind, C, lbd) {
  K <- nrow(C)
  cvec <- c(1, rep(0, 2*K))
  Amat <- -cvec
  Amat <- rbind(Amat, c(-1, rep(1, 2*K)))
  tmp_constr <- C %x% t(c(1,-1))
  Amat <- rbind(Amat, cbind(-1 * lbd, rbind(tmp_constr, -tmp_constr)))
  tmp_vec <- rep(0, K)
  tmp_vec[col_ind] <- 1
  bvec <- c(0, 0, tmp_vec, -tmp_vec)

  lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  while (length(lpResult) == 0) {
    cat("The penalty lambda =", lbd, "is too small and increased by 0.01...\n")
    lbd <- lbd + 0.01
    Amat[-(1:2),1] <- lbd
    lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2*K, 2)
  return(lpResult[ind] - lpResult[ind + 1])
}

#' This function estimate omega
#'
#' For a given lbd and C, solve the C^{-1}.
#' Require: C should be symmetric and square.
#' @param lbd list of lambda values
#' @param C K by K estimated C_hat
#' @return omega
#' @export
EstOmega <- function(lbd, C) {
  K <- nrow(C)
  omega <- matrix(0, K, K)
  for (i in 1:K) {
    omega[,i] <- Solve_row(i, C, lbd)
  }
  return(omega)
}


######################################################################
#####                                                            #####
##### This is the code for estimating the pure node set          #####
#####                                                            #####
######################################################################

#' Calculate the maximal absolute value for each row of the given matrix.
#'
#' @param  Sigma p by p matrix
#' @return length p vector
#' @export
FindRowMax <- function(Sigma) {

  M <- rep(0, nrow(Sigma))
  for (i in 1:nrow(Sigma)) {
    Mi <- max(abs(Sigma[i, -i]))   # remove the diagonal element
    M[i] <- Mi
  }
  return(M)
}

#' Calculate indices of ith row such that the absolute values of these indices
#' are within 2*delta from the maximal absolute value {@code M} of this row.
#'
#' @param i integer denoting for the ith row.
#' @param M the maximal absolute value of the ith row.
#' @param vector the entries of this row.
#' @param delta numerical constant.
#' @return a vector of indices.
#' @export
FindRowMaxInd <- function(i, M, vector, delta) {
  lbd <- 2 * delta
  indices <- which(M - abs(vector) <= lbd)
  return(indices[which(indices != i)])
}


#' For given row, check if it is a pure node by interactively checking the nodes in Si. Return TRUE if the given row corresponds to a pure variable.
#'
#' @param Sigma p by p matrix.
#' @param rowInd integer index
#' @param Si vector of indices.
#' @param Ms vector of largest absolute values of each rows in Si
#' @param delta numerical constant.
#' @return TRUE or FALSE.
#' @export
TestPure <- function(Sigma, rowInd, Si, Ms, delta) {


  for (i in 1:length(Si)) {
    j <- Si[i]
    if (abs(abs(Sigma[rowInd, j]) - Ms[j]) > 2 * delta)
      return(FALSE)
  }
  return(TRUE)
}

#' Union all the groups which have common variables
#'
#' @param groupList groups of clustering
#' @param Union flag for using union or intersection for p sets of pure nodes.
#' @return TRUE or FALSE.
#' @export
#'
UnionGroup <- function(groupList, Union) {
  newGroup <- list(groupList[[1]])
  for (i in 2:length(groupList)) {
    groupi <- groupList[[i]]
    inGroupFlag <- FALSE
    for (j in 1:length(newGroup)) {
      groupj <- newGroup[[j]]
      if (length(intersect(groupi, groupj)) != 0) {
        inGroupFlag <- TRUE
        if (Union) {
          newGroup[[j]] <- union(groupi, groupj)
        } else {
          newGroup[[j]] <- intersect(groupi, groupj)
        }
        break
      }
    }
    if (!inGroupFlag)
      newGroup <- append(newGroup, list(groupi))
  }
  return(newGroup)
}

#' Check if all the sets are disjoint
#'
#' @param groupList groups of clustering
#' @return TRUE or FALSE.
#' @export
#'
CheckDisjoint <- function(groupList) {
  if (length(groupList) == 1)
    return(TRUE)
  for (i in 1:(length(groupList) - 1)) {
    groupi <- groupList[[i]]
    for (j in (i+1):length(groupList)) {
      groupj <- groupList[[j]]
      if (length(intersect(groupi, groupj)) != 0)
        return(FALSE)
    }
  }
  return(TRUE)
}

#' Estimate list of pure node indices for given {@code Sigma} and {@delta}.
#'
#' @param Sigma p by p sample covariance matrix.
#' @param delta numerical constant.
#' @param Ms the largest absolute values of each row of Sigma.
#' @param Union flag for using union or intersection for p sets of pure nodes.
#' @return a list including:
#' \itemize{
#' \item the list of the estimated pure node indices
#' \item the vector of the estimated pure node indices
#' }
#' @export
#'
FindPureNode = function(Sigma, delta, Ms, Union) {
  G <- list()
  vec <- c()
  # Calculate the Si for those rows which pass the checking
  for (i in 1:nrow(Sigma)) {
    Si <- FindRowMaxInd(i, Ms[i], Sigma[i,], delta)
    if (length(Si) == 0)
      next
    pureFlag <- TestPure(Sigma, i, Si, Ms, delta)
    if (pureFlag) {
      Gi <- c(Si, i)
      G <- append(G, list(Gi))
      names(G)[length(G)] <- i
    }
  }
  while(!CheckDisjoint(G))
    G <- UnionGroup(G, Union)
  return(list(pureInd = G, pureVec = unlist(G)))
}



#' Estimate the sign subpartition of pure node sets. If there is an element of a list is empty, then a empty list will be put in that position
#'
#' @param pureList list of pure node indices (Example: list(c(1,2,3),c(4,5,6,7)))
#' @param Sigma p by p sample covariance
#' @return list of sign subpartition of pure node indices.
#' Example: list(list(c(1,2),3),list(c(4,7),c(5,6)))
#' @export
#'
FindSignPureNode <- function(pureList, Sigma) {
  signPureList <- list()
  for (i in 1:length(pureList)) {
    purei <- pureList[[i]]
    if (length(purei) == 1) {
      signPureList[[i]] <- list(pos = purei[1], neg = list())
      next
    }
    firstPure <- purei[1]
    pos <- firstPure
    neg <- c()
    for (j in 2:length(purei)) {
      purej <- purei[j]
      if (Sigma[firstPure, purej] < 0)
        neg <- c(neg, purej)
      else
        pos <- c(pos, purej)
    }
    if (length(pos) == 0) {
      pos <- list()
    } else {
      pos <- pos
    }
    if (length(neg) == 0) {
      neg <- list()
    } else {
      neg <- neg
    }
    signPureList[[i]] <- list(pos = pos, neg = neg)
  }
  return(signPureList)
}


#' Recover the estimated submatrix A_I by given the pure node group.
#'
#' @param estGroupList list of group indices of the pure node with sign.
#' @param p number of rows
#' @return p by K matrix.
#' @export
#'
RecoverAI <- function(estGroupList, p) {

  K <- length(estGroupList)
  A <- matrix(0, p, K)
  for (i in 1:K) {
    groupi <- estGroupList[[i]]
    groupi1 <- groupi[[1]]       # group with positive sign
    for (j in groupi1)
      A[j,i] <- 1
    groupi2 <- groupi[[2]]
    if (length(groupi2) != 0) {
      for (k in groupi2) {
        A[k,i] <- -1
      }
    }
  }
  return(A)
}

#' Use the given {@code delta} to calculate the fitted AI, pure variable indices in list form and vector form. Also return estimated Y and C for the following Dantzig estimation.
#'
#' @param Sigma p by p covariance matrix.
#' @param delta delta to be used.
#' @param Union flag for using union or intersection for p sets of pure nodes.
#' @return the list including:
#' \itemize{
#' \item AI: the p by K estimated AI
#' \item pureVec: vector of the indices of estimated pure variables
#' \item pureSignInd: list of the indices of estimated pure variables
#' }
#' @export
#'
EstAI <- function(Sigma,delta,Union) {

  Ms <- FindRowMax(Sigma)
  resultPure <- FindPureNode(Sigma, delta, Ms, Union)
  estPureIndices <- resultPure$pureInd
  estPureVec <- resultPure$pureVec
  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, nrow(Sigma))
  return(list(AI = AI, pureVec = estPureVec, pureSignInd = estSignPureIndices))
}



################################################################################
######                                                                   #######
######   These are the code which convert the original problem into LP   #######
######                                                                   #######
################################################################################

#'For a given ith row of the matrix Omega and {@code vec}, find the vec_constr
#'
#'
#'such that t(vec_constr) %*% par_vec = 't(vec) %*% (ci1(+),ci1(-),ci2(+),ci2(-),ci3(+),ci3(-)).
#'
#' @param i is an index of each row, 1 <= i <= n.
#' @param vec length 2*p vector containing the values to be modified.
#' @return length p*(p+1) vector.
#'  eg. |C21| + |C22| + |C23| = t(c(1,1,1,1,1,1) %*% c(c21(+),c21(-),c22(+),c22(-),c23(+),c23(-))
#'  = t(c(0,0,1,-1,0,0,1,-1,1,-1,0,0)) %*% par_vec
#' @export
#'
setConstr <- function(i, vec) {
  # Returns:
  #   length p*(p+1) vector.
  #   eg: |C21| + |C22| + |C23| = t(c(1,1,1,1,1,1) %*% c(c21(+),c21(-),c22(+),c22(-),c23(+),c23(-))
  #                           = t(c(0,0,1,-1,0,0,1,-1,1,-1,0,0)) %*% par_vec
  K <- length(vec) / 2
  constr <- c()
  for (j in 1:K) {
    if (j <= i - 1) {
      tmpConstr <- rep(0, 2 * (K - j + 1))
      tmpInd <- i - j + 1
      tmpConstr[(2 * tmpInd - 1):(2 * tmpInd)] = vec[(2 * j - 1):(2 * j)]
      constr <- c(constr, tmpConstr)
    } else {
      tmpConstr <- vec[(2 * j - 1):(2 * K)]
      constr <- c(constr, tmpConstr)
      break
    }
  }
  constr <- c(constr, rep(0, K * (K + 1) - length(constr)))
  return(constr)
}


#' Get the Amat constraint matrix of Inf-1 norm, i.e. the l_1 norm of each row.
#'
#' @param p numerical constant denoting the dimension of Omega.
#' @return p by p(p+1) matrix.
#' @export
#'
lInfConstr <- function(p) {
  constrMat <- matrix(0, p, p * (p + 1))
  vec <- rep(1, 2 * p)
  for (i in 1:p) {
    constrMat[i, ] = setConstr(i, vec)
  }
  return(constrMat)
}


#' Get the Amat constraint matrix for the product Omega %*% C based on
#' diagonal and off-diagonal such that:
#'       t(A_i) %*% par_vec <- t(Omega_i) %*% C_i,
#' for i belonging to the diagonal set. Similar for the off-diagonal,
#' the order for each row of off-diagonal set is from the first row to
#' the last row and ordered by row.
#'
#' @param C K by K matrix.
#' @return a list including:
#' \itemize{
#' \item the constraint matrices of diagonal and offdiagonal.
#' }
#' @export
#'
prodConstr <- function(C) {

  K <- nrow(C)
  constrDiag <- matrix(0, K, K * (K + 1))
  constrOffDiag <- matrix(0, K * (K - 1) / 2, K * (K + 1))
  counter <- 0
  for (i in 1:K) {
    for (j in i:K) {
      vec <- rep(C[ ,j], each = 2) * c(1,-1)
      if (j == i) {
        constrDiag[j, ] <- setConstr(j, vec)
      } else {
        counter <- counter + 1
        constrOffDiag[counter, ] <- setConstr(i, vec)
      }
    }
  }
  return(list(diag = constrDiag, off = constrOffDiag))
}


#' Obtain the Amat and bvec for a given C and lbd.
#'
#' @param lbd lamdba value
#' @param C K by K estimated C_hat
#' @return a list including:
#' \itemize{
#' \item Amat
#' \item bvec
#' }
#' @export
genConstr <- function(lbd, C) {
  K <- nrow(C)
  A <- diag(-1, nrow = 1 + K * (K + 1))
  constr1 <- cbind(-1, lInfConstr(K))
  A <- rbind(A, constr1)
  b <- c(rep(0, K * (K + 1) + 1), rep(0, nrow(constr1)))
  result <- prodConstr(C)
  constr2 <- rbind(cbind(-lbd, result$diag), cbind(-lbd, -result$diag))
  b <- c(b, rep(c(1,-1), each = K))
  constr3 <- rbind(cbind(-lbd, result$off), cbind(-lbd, -result$off))
  b <- c(b, rep(0, 2 * nrow(result$off)))
  A <- rbind(A, constr2, constr3)
  return(list(A = A, b = b))
}



#' Use LP to solve problem:
#'
#' @param y K by 1 vector.
#' @param lbd positive contant.
#' @return K by 1 vector
#' @export
LP <- function(y, lbd) {

  # beta^+ - beta^- \leq lbd + y
  # - beta^+ + beta^- \leq lbd - y
  # beta^+ \ge 0; beta^- \ge 0

  K <- length(y)
  cvec <- rep(1, 2 * K)
  bvec <- c(lbd + y, lbd - y, rep(0, 2 * K))
  C <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    indices <- c(i,i + K)
    C[i,indices] = c(1,-1)
  }
  Amat <- rbind(C, -C, diag(-1, nrow=2 * K))
  LPsol <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)]
  return(beta)
}



################################################################################
######                                                                   #######
######   These are the code for Parameter tuning and Cross validation    #######
######                                                                   #######
################################################################################


#' Cross validation for choosing delta. For each delta from the given grids, first split the data into two datasets.
#' Obtain I, AI and C. Then calculate AICAI^T and choose delta which minimizes the criterion: AICAI^t - Sigma(dataset 2)
#'
#' @param X n by p matrix.
#' @param deltaGrids vector of numerical constants.
#' @param Union TRUE or FALSE for estimating the pure variables. (For details, see doc of FindPureNode.)
#' @param diagonal TRUE or FALSE for the diagonal structure of C.
#' @return the selected optimal delta
#' @export
CV_Delta <- function(X, deltaGrids, Union, diagonal) {

  sampInd <- sample(dim(X)[1], floor(dim(X)[1] / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- t(X1) %*% X1 / nrow(X1)
  Sigma2 <- t(X2) %*% X2 / nrow(X2)
  Ms <- FindRowMax(Sigma1)
  loss <- c()
  for (i in 1:length(deltaGrids)) {
    resultFitted <- CalFittedSigma(Sigma1, deltaGrids[i], Ms, Union, diagonal)
    fittedValue <- resultFitted$fitted
    estPureVec <- resultFitted$pureVec
    if (is.null(dim(fittedValue)) && fittedValue == -1) {
      loss[i] <- Inf
    } else {
      subSigma2 <- Sigma2[estPureVec, estPureVec]
      denom <- length(estPureVec) * (length(estPureVec)-1)
      loss[i] <- 2 * offSum(subSigma2, fittedValue) / denom
    }
  }
  return(deltaGrids[stats::median(which(loss == min(loss)))])
}


#' Calculate the fitted value of A_ICA_I^T for given Sigma and delta.
#'
#' @param Sigma p by p covariance matrix.
#' @param delta given parameter.
#' @param Ms the calculated maximal values of Sigma per row.
#' @param Union TRUE or FALSE.
#' @param diagonal TRUE or FALSE.
#' @return a list including:
#' \itemize{
#' \item pureVec: vector of the indices of estimated pure variables.
#' \item fitted: fitted value of A_ICA_I^T
#' }
#
#' @export
CalFittedSigma <- function(Sigma, delta, Ms, Union, diagonal) {

  resultPureNode <- FindPureNode(Sigma, delta, Ms, Union)
  estPureIndices <- resultPureNode$pureInd
  if (singleton(estPureIndices))
    return(list(pureVec = NULL, fitted = -1))
  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, nrow(Sigma))
  C <- EstC(Sigma, AI, diagonal)
  if (length(estPureIndices) == 1) {
    fitted <- -1
  } else {
    subAI <- AI[resultPureNode$pureVec, ]
    fitted <- subAI %*% C %*% t(subAI)
  }
  return(list(pureVec = resultPureNode$pureVec, fitted = fitted))
}

#' Use Cross-validation to select lambda for estimating Omega.
#'
#' Split the data into two parts. Estimating C on two datasets. Then, for each lambda, calculate Omega on the first dataset and calculate the loss on the second dataset.
#' Find the lambda which gives the smallest loss of <C,Omega> - log(det(Omega))
#' @param X n by p matrix.
#' @param lbdGrids vector of numerical constants.
#' @param AI p by K matrix
#' @param pureVec vector of indices of pure variables.
#' @param diagonal TRUE or FALSE for the diagonal structure of C.
#' @return the selected optimal lambda
#
#' @export
CV_lbd <- function(X, lbdGrids, AI, pureVec, diagonal) {

  sampInd <- sample(nrow(X), floor(nrow(X) / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  print("calculating Sigma1")
  Sigma1 <- t(X1) %*% X1 / nrow(X1)
  print("calculating Sigma2")
  Sigma2 <- t(X2) %*% X2 / nrow(X2)
  print("calculating C1")
  C1 <- EstC(Sigma1, AI, diagonal)
  print("calculating C2")
  C2 <- EstC(Sigma2, AI, diagonal)
  loss <- c()
  for (i in 1:length(lbdGrids)) {
    Omega <- EstOmega(lbdGrids[i], C1)
    det_Omega <- det(Omega)
    if (det_Omega <= 0) {
      loss[i] <- Inf
    } else {
      loss[i] <- sum(Omega * C2) - log(det_Omega)
    }
    print(paste0('lambda: ',lbdGrids[i]," det_Omega: ",det_Omega," loss ", loss[i]))
  }

  return(lbdGrids[which.min(loss)])
}



#' Check if an element is in a list. If it does exist in some group, return the group
#' index in that list and its sublist index. Otherwise, return c(0,0).
#'
#' @param element a variable
#' @param groupList list of variable indices for each group
#' @return index of element in the list and its sublist, If the element is not in the list, return c(0,0)
#
#' @export
checkElement <- function(element, groupList) {

  for (i in 1:length(groupList)) {
    for (j in 1:length(groupList[[i]])) {
      if (element %in% groupList[[i]][[j]])
        return(c(i,j))
    }
  }
  return(c(0,0))
}


#' For a given integer K, generate its sign permutation
#'
#' @param n an integer K
#' @return a 2^K by K matrix containing all sign permutations
#' @export
signPerm <- function(n) {
  signPerms <- rep(1, n)
  for (i in 1:n) {
    allCombs <- gtools::combinations(n, i, 1:n)
    for (j in 1:nrow(allCombs)) {
      thisPerm <- rep(1, n)
      thisPerm[allCombs[j, ]] <- -1
      signPerms <- rbind(signPerms, thisPerm)
    }
  }
  return(signPerms)
}


#' For given matrix A, find the row indices which correspond to pure nodes
#'
#' @param A p by K matrix.
#' @return vector of pure row indices
#' @export
pureRowInd <- function(A) {

  pureVec <- c()
  for (i in 1:ncol(A)) {
    pureVec <- c(pureVec, which(abs(A[ ,i]) == 1))
  }
  return(pureVec)
}



#' For the given pure node set, find the best column permutation matrix
#' Aperm of A such that |Aperm - B|_F is minimal for the target B
#'
#' @param A generic matrix.
#' @param B generic matrix with the same dimension of A.
#' @return list of column permutation and sign permutation
#' @export
getOptPerm <- function(A, B) {
  K <- ncol(A)
  allPerm <- gtools::permutations(K, K, 1:K)
  signPerms <- signPerm(K)
  prevLoss <- norm(A - B, "F")
  optPerm <- vector("list", length <- 2)
  for (i in 1:nrow(allPerm)) {
    # skip the identity permutation
    permi <- methods::as(as.integer(allPerm[i, ]), "pMatrix")
    permutatedA <- A %*% permi
    for (j in 1:nrow(signPerms)) {
      signPermj <- signPerms[j, ]
      # perform sign permutation along columns of newA
      newA <- t(t(permutatedA) * signPermj)
      currLoss <- norm(newA - B, "F")
      if (currLoss <= prevLoss) {
        optPerm[[1]] <- permi
        optPerm[[2]] <- signPermj
        prevLoss <- currLoss
      }
    }
  }
  return(optPerm)
}


#' Find the optimal permutated matrix of A such that |Aperm-B|_F is minimized.
#'
#' @param  A,B two matrices with the same dimension.
#' @return Permutated A.
#' @export
permA <- function(A, B) {
  if (sum(dim(A) == dim(B)) == 2) {
    pureInd <- pureRowInd(B)
    optPerm <- getOptPerm(A[pureInd, ], B[pureInd, ])
    permutatedA <- A %*% optPerm[[1]]
    A <- t(t(permutatedA) * optPerm[[2]])
  }
  return(A)
}



#' Calculate specificity and sensitivity based on group partition
#'
#' @param  p the total number of variables
#' @param estGroup list of variable indices for estimated group
#' @param trueGroup list of variable indices for true group
#' @return Specificity and Sensitivity
#' @export
calSPandSN <- function(p, estGroup, trueGroup) {

  TP <- TN <- FP <- FN <- TPS <- TNS <- FNS <- FPS <- 0
  estTable <- trueTable <- matrix(0, nrow=p, ncol=2)
  for (k in 1:p) {
    estTable[k,] <- checkElement(k,estGroup)
    trueTable[k,] <- checkElement(k,trueGroup)
  }
  for (i in 1:(p-1)) {
    flagiEst <- estTable[i,]
    flagiTrue <- trueTable[i,]
    for (j in (i+1):p) {
      flagjEst <- estTable[j,]
      flagjTrue <- trueTable[j,]
      if (flagiTrue[1] * flagjTrue[1] > 0 && flagjTrue[1] == flagiTrue[1]) {
        if (flagiEst[1] * flagjEst[1] > 0 && flagiEst[1] == flagjEst[1]) {
          TP <- TP + 1
          if (flagiTrue[2] == flagjTrue[2]) {
            if (flagiEst[2] == flagjEst[2]) {
              TPS <- TPS + 1
            } else {
              FNS <- FNS + 1
            }
          } else {
            if (flagiEst[2] == flagjEst[2]) {
              FPS <- FPS + 1
            } else {
              TNS <- TNS + 1
            }
          }
        } else {
          FN <- FN + 1
        }
      } else {
        if (flagiEst[1] * flagjEst[1] > 0 && flagiEst[1] == flagjEst[1]) {
          FP <- FP + 1
        } else {
          TN <- TN + 1
        }
      }
    }
  }
  ifelse (TN + FP == 0, SP <- 0, SP <- TN / (TN + FP))
  ifelse (TP + FN == 0, SN <- 0, SN <- TP / (TP + FN))
  ifelse (FNS + TPS == 0, SNS <- 0, SNS <- TPS / (FNS + TPS))
  ifelse (TNS+FPS == 0, SPS <- 0, SPS <- TNS / (FPS + TNS))
  return(c(SP = SP, SN = SN, SPS = SPS, SNS = SNS))
}



#' Recover group structure based given p by K matrix A and perform thresholding
#'
#' @param  A estimated matrix. thresh constant > 0.
#' @return list of group indices with sign subpartition
#' @export
recoverGroup <- function(A) {

  Group <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    Group[[i]] <- list(pos = posInd, neg = negInd)
  }
  return(Group)
}


#' Calculate the FNR and FPR for given p by K matrix A.
#'
#' @param  estA,trueA two matrices having same dims.
#' @return a list including:
#' \itemize{
#' \item list of error: FNR and FPR.
#' \item errorInd: the row indices where errors occur.
#' }
#' @export
calFNRandFPR <- function(estA,trueA) {

  if (sum(dim(estA) == dim(trueA)) != 2) {
    # If estA and A don't have the same dimension, we return -1.
    #    cat("Two matrices don't have the same dimension!\n")
    return(list(error    = c(FPR = -1, FNR = -1, FPSR=-1, FNSR=-1),
                errorInd = c()))
  }
  errorInd <- c()
  p <- nrow(estA)
  K <- ncol(estA)
  TP <- FP <- TN <- FN <- FPS <- TPS <- TNS <- FNS <- 0
  for (i in 1:p) {
    for (j in 1:K) {
      if (trueA[i,j] != 0) {
        TN <- TN + 1
        if (estA[i,j] == 0) {
          errorInd <- c(errorInd, i) # record the row indices where errors occur
          FN <- FN + 1
        } else {  ### calculate the sign error rate
          if (trueA[i,j] > 0) {
            TPS <- TPS + 1
            if (estA[i,j] < 0)
              FNS <- FNS + 1
          } else {
            TNS <- TNS + 1
            if (estA[i,j] > 0)
              FPS <- FPS + 1
          }
        }
      } else {
        TP <- TP + 1
        if (estA[i,j] != 0) {
          FP <- FP + 1
          errorInd <- c(errorInd, i)
        }
      }
    }
  }
  if (TPS == 0)
    FNSR <- 0
  if (TNS == 0)
    FPSR <- 0
  return(
    list(error    = c(FPR = FP / TP, FNR = FN / TN,
                      FPSR = FPS / TNS, FNSR = FNS / TPS),
         errorInd = unique(errorInd))
  )
}



#' Calculate the matrix sup norm, frobenius norm(scaled) and L_inf norm.
#'
#' @param  A matrix
#' @return vector of three numerical results of c(frob, sup, L_inf).
#' @export
calArates <- function(A) {

  frob <- norm(A, "F") / sqrt(nrow(A))
  max <- norm(A, "M")
  inf <- sum(abs(A)) / nrow(A)
  return(c(sup = max, l1 = inf, l2 = frob))
}


#' Check if there exists an element of the given list has length equal to 1
#'
#'
#' If exists at least one, return TRUE; otherwise return FALSE
#'
#' @param  estPureIndices list of pure indices
#' @return True or False
#' @export
singleton <- function(estPureIndices) {

  for (i in 1:length(estPureIndices)) {
    if (length(estPureIndices[[i]]) == 1)
      return(T)
  }
  return(F)
}


#' Infer error between estimated and true matrices
#'

#' @param estA estimated matrix
#' @param trueA true matrix
#' @param flagPerm permutation or not
#' @return a list of (specificity and sensitivity, error,matrix sup norm, frobenius norm(scaled) and L_inf norm)
#' @export
getErrors <- function(estA, trueA, flagPerm = TRUE) {
  estGroup <- recoverGroup(estA)
  trueGroup <- recoverGroup(trueA)
  speAndSen <- calSPandSN(nrow(trueA), estGroup, trueGroup)
  if (flagPerm)
    estAperm <- permA(estA, trueA)
  else
    estAperm <- estA
  fr <- calFNRandFPR(estAperm, trueA)$error
  if (fr[1] == -1) {
    rates <- rep(-1, 3)
  } else {
    rates <- calArates(estAperm - trueA)
  }
  return(list(SPandSN = speAndSen, FNRandFPR = fr,rates = rates))
}



#' Normalize matrix based on threshold
#'
#' Threshold the estimated {@code A} based on the given {@code mu}. If {@code scale} is true,
#' then normalize each row of A such that the l-1 norm of each row is not larger than 1.
#'
#' @param A estimated matrix
#' @param mu threshold
#' @param scale normalize or not
#' @return normalized matrix
#' @export
threshA <- function(A, mu, scale = FALSE) {
  if (scale) {
    scaledA <- A /rowSums(abs(A))
  } else {
    scaledA <- A
  }
  for (i in 1:nrow(A)) {
    colInd <- abs(scaledA[i, ]) <= mu
    scaledA[i,colInd] = 0
  }
  return(scaledA)
}


#' Calculate the sum of squares of the upper off-diagonal elements of two matrices (require: M and N have the same dimensions)
#'
#' @param M Matrix M
#' @param N Matrix N
#' @return sum of squares of the upper off-diagonal elements of two matrices
#' @export
offSum <- function(M, N) {
  tmp <- M-N
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}
