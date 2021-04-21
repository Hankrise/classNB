#' Parameter initialisation 
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of mean (optional).
#' @param var a vector of variance (optional).
#' @param initPr a vector of initial probability of HMM (optional).
#' @param transPr a matrix of transition probability matrix of HMM (optional).
#' @return A list of 5 objets.
#' \describe{
#' \item{\code{esAvg}}{a vector of the initialised parameter of mean.}
#' \item{\code{esVar}}{a vector of the initialised parameter of variance.}
#' \item{\code{eSize}}{a vector of the initialised parameter of size.}
#' \item{\code{eProb}}{a vector of the initialised parameter of probability.}
#' \item{\code{initPr}}{a vector of the initialised parameter of initial probability of HMM.}
#' \item{\code{transPr}}{a matrix of the initialised parameter of transition probability matrix of HMM.}
#' }
HMMNB.init <- function(X, nbS,avg, var, initPr, transPr){
  if(is.null(initPr)|is.null(transPr)){
    initPr = rep(1/nbS, nbS)
    transPr = 0.5*diag(nbS) + array(0.5/(nbS),c(nbS,nbS))
  }

  if(is.null(avg)|is.null(var)){
    esAvg=1;esVar=0;
    while(any(esAvg>esVar)){
      clust <- kmeans(x = as.vector(t(X)), centers = nbS)
      esAvg = as.vector(clust$centers)
      ordAvg=order(esAvg)
      esAvg = esAvg[ordAvg]
      esVar = clust$withinss / (clust$size - 1)
      esVar = esVar[ordAvg]
    }
  } else {
    esAvg=avg
    esVar=var
  }
  eSize = esAvg^2/(esVar-esAvg)
  eProb = esAvg/esVar

  return(list(esAvg=esAvg,esVar=esVar,eSize=eSize,eProb=eProb,initPr=initPr,transPr=transPr))
}
#' Perform inference of HMM negative binomial models
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of initial value of mean (optional).
#' @param var a vector of initial value of variance (optional).
#' @param initPr a vector of initial probability of HMM (optional).
#' @param transPr a matrix of transition probability matrix of HMM (optional).
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_CF algorithm (optional).
#' @param EPSILON a value for the threshold used for the stopping criteria for the mixNB_CF algorithm (optional).
#' @param algo a string of specifying the approach of M-step in EM algorithm, available choices are "EM-FL", "GEM-NM", "GEM-BFGS", "GEM-DFP".
#' @param accelerate a string of specifying the acceleration method of EM algorithm, available choices are NULL(no accleration), "SQEM", "DAAREM".
#' @return A list of 3 objets.
#' \describe{
#' \item{\code{parameter}}{a list with the inferred parameters: Size, Probability, mean, variance, initial probability and transition probability of HMM.}
#' \item{\code{emisPr}}{a matrix with emission probability associated to each states in column and each position in row.}
#' \item{\code{postPr}}{a matrix with posterior probability associated to each series in column and each position in row.}
#' \item{\code{iter}}{a number of iterations while algorithm stop.}
#' }
#' @export
#' @examples
#' data(toyexample)
#' # Perform inference of HMM negative binomial models
#' resHMMNB <- HMMNB(X = toydata, nbS = 3,algo="EM-FL")
#' @seealso \code{\link{HMMNB_NS}}, \code{\link{HMMNB_NS}}
HMMNB <- function(X, nbS, avg=NULL, var=NULL, initPr=NULL, transPr=NULL, iterMax=15000, EPSILON=1e-7,algo="EM-FL",accelerate=NULL){
  init=HMMNB.init(X, nbS,avg, var, initPr, transPr)

  nbT=length(X)
  esAvg=init$esAvg
  esVar=init$esVar
  eSize=init$eSize
  eProb=init$eProb
  initPr=init$initPr
  transPr=init$transPr
  ones_nbS=as.matrix(rep(1,nbS))
  ones_nbT=as.matrix(rep(1,nbT))
  ones_nbTS=ones_nbT%*%t(ones_nbS)
  ones_nbST=ones_nbS%*%t(ones_nbT)
  ones_nbSS=ones_nbS%*%t(ones_nbS)

  if(is.null(accelerate)){
    #iteration starts
    for(iter in 1:iterMax){
      prevParam = c(esAvg, eSize, initPr, c(transPr)) # define iterative parameters in an arrray
      Nmat = eSize%*%t(ones_nbT) # matrix of size
      Mmat = esAvg%*%t(ones_nbT) # matrix of mean
      Xmat = ones_nbS%*%t(X)     # matrix of sample X
      ## Estep --------
      Pmat = dnbinom(Xmat,size=Nmat,mu=Mmat) # matrix of emission Probability

      Fpr = matrix(NA, nbS, nbT)
      Apr = rep(NA, nbT)
      Apr[1]  = sum(initPr*Pmat[,1])
      Fpr[,1]  = initPr*Pmat[,1]/Apr[1]
      for(tt in 2:nbT){
        Fpr.tmp = rep(0,nbS)
        for(l in 1:nbS){
          Fpr.tmp[l] = Pmat[l,tt]*sum(Fpr[,tt-1]*transPr[,l])
        }
        Apr[tt]  = sum(Fpr.tmp)
        Fpr[,tt] = Fpr.tmp/Apr[tt]
      }
      resFwd = list(Fpr=Fpr, logLik=sum(log(Apr)))
      Fpr = resFwd$Fpr
      logLik = resFwd$logLik
      nbT = ncol(Fpr); nbS = nrow(Fpr) # backward
      postPr = matrix(NA, nbS, nbT)
      Gpr = matrix(NA, nbS, nbT)
      postPr[,nbT] = Fpr[,nbT]
      for(tt in (nbT-1):1){
        Gpr[,tt+1]  = Fpr[,tt] %*% transPr
        postPr[,tt] = Fpr[,tt] *(transPr %*% (postPr[,tt+1]/Gpr[,tt+1]))#sapply(1:nbS, function(l) sum(trans[l,] * postPr[,tt+1] / Gpr[,tt+1]))#
      }
      resBcd  = list(Gpr = Gpr, postPr = postPr)
      postPr  = resBcd$postPr
      Gpr = resBcd$Gpr

      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
        esAvg = sort(eSize*(1-eProb)/eProb);  # new value of mean
        esVar = eSize*(1-eProb)/eProb^2;  # new value of variance
      }else if(algo=="GEM-NM"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          eParam = Nelder.Mead(objFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-BFGS"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = BFGS(objFunc, devFunc, par)$xmin
          # eParam = optim(par,objFunc,devFunc,method = "BFGS")$par
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-DFP"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = DFP(objFunc, devFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }

      Pi.tmp = matrix(NA, nbS, nbS)
      eParam = NULL
      for(k in 1:nbS){
        for(l in 1:nbS){
          Pi.tmp[k,l] = transPr[k,l]*sum( Fpr[k,1:(nbT-1)]*postPr[l,2:nbT]/Gpr[l,2:nbT] )
        }
        Mean.tmp = sum(postPr[k,] * X) / sum(postPr[k,])
        Var.tmp  = sum(postPr[k,] * (X - Mean.tmp)^2) / sum(postPr[k,])
        MV.tmp = c(Mean.tmp, Var.tmp)
        eParam = rbind(eParam, MV.tmp)
      }
      transPr = Pi.tmp/rowSums(Pi.tmp)
      initPr = postPr[,1]
      newParam = c(esAvg, eSize, initPr, c(transPr)) # update iterative parameters in an arrray
      # throw warning when any element in the iterative parameters can't be calculated
      if(any(is.na(esAvg)) || any(is.na(esAvg))){
        warning("HMM EM - NA is appeared!")
        esAvg=prevParam[1:nbS]; eSize=prevParam[(nbS+1):(2*nbS)];
        break
      }
      D = sum((newParam - prevParam)^2) # when precision reach the desired level, stop iteration, otherwise reach the maximum limit
      if(D < EPSILON){break}
    }
    cat("iterStop = ", iter, "\n")
    # return with a list
    parameter=list(); parameter[["eSize"]]=eSize; parameter[["eProb"]]=eProb; parameter[["esAvg"]]=esAvg;
    parameter[["esVar"]]=esVar; parameter[["initPr"]]=as.numeric(initPr); parameter[["transPr"]]=as.numeric(transPr)
    structure(list(parameter=parameter, emisPr=Pmat, postPr=postPr,iter=iter), class="HMMNB")

  }else if(accelerate=="SQEM"){
    step=1
    #iteration starts
    for(iter in 1:iterMax){
      prevParam = c(esAvg, eSize, initPr, c(transPr)) # define iterative parameters in an arrray

      ## 1. theta1 ---------------------------------------------------
      Nmat = eSize%*%t(ones_nbT) # matrix of size
      Mmat = esAvg%*%t(ones_nbT) # matrix of mean
      Xmat = ones_nbS%*%t(X)     # matrix of sample X
      ## Estep --------
      Pmat = dnbinom(Xmat,size=Nmat,mu=Mmat) # matrix of emission Probability
      nbT = ncol(Pmat); nbS = nrow(Pmat) # forward
      Fpr = matrix(NA, nbS, nbT)
      Apr = rep(NA, nbT)
      Apr[1]  = sum(initPr*Pmat[,1])
      Fpr[,1]  = initPr*Pmat[,1]/Apr[1]
      for(tt in 2:nbT){
        Fpr.tmp = rep(0,nbS)
        for(l in 1:nbS){
          Fpr.tmp[l] = Pmat[l,tt]*sum(Fpr[,tt-1]*transPr[,l])
        }
        Apr[tt]  = sum(Fpr.tmp)
        Fpr[,tt] = Fpr.tmp/Apr[tt]
      }
      resFwd = list(Fpr=Fpr, logLik=sum(log(Apr)))
      Fpr = resFwd$Fpr
      logLik = resFwd$logLik
      nbT = ncol(Fpr); nbS = nrow(Fpr) # backward
      postPr = matrix(NA, nbS, nbT)
      Gpr = matrix(NA, nbS, nbT)
      postPr[,nbT] = Fpr[,nbT]
      for(tt in (nbT-1):1){
        Gpr[,tt+1]  = Fpr[,tt] %*% transPr
        postPr[,tt] = Fpr[,tt] *(transPr %*% (postPr[,tt+1]/Gpr[,tt+1]))#sapply(1:nbS, function(l) sum(trans[l,] * postPr[,tt+1] / Gpr[,tt+1]))#
      }
      resBcd  = list(Gpr = Gpr, postPr = postPr)
      postPr  = resBcd$postPr
      Gpr = resBcd$Gpr
      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
        esAvg = sort(eSize*(1-eProb)/eProb);  # new value of mean
        esVar = eSize*(1-eProb)/eProb^2;  # new value of variance
      }else if(algo=="GEM-NM"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          eParam = Nelder.Mead(objFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-BFGS"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = BFGS(objFunc, devFunc, par)$xmin
          # eParam = optim(par,objFunc,devFunc,method = "BFGS")$par
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-DFP"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = DFP(objFunc, devFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }
      Pi.tmp = matrix(NA, nbS, nbS)
      eParam = NULL
      for(k in 1:nbS){
        for(l in 1:nbS){
          Pi.tmp[k,l] = transPr[k,l]*sum( Fpr[k,1:(nbT-1)]*postPr[l,2:nbT]/Gpr[l,2:nbT] )
        }
        Mean.tmp = sum(postPr[k,] * X) / sum(postPr[k,])
        Var.tmp  = sum(postPr[k,] * (X - Mean.tmp)^2) / sum(postPr[k,])
        MV.tmp = c(Mean.tmp, Var.tmp)
        eParam = rbind(eParam, MV.tmp)
      }
      transPr = Pi.tmp/rowSums(Pi.tmp)
      initPr = postPr[,1]
      esAvg1 = esAvg; eSize1 = eSize; initPr1 = initPr; transPr1 = transPr
      newParam1 = c(esAvg1, eSize1, initPr1, c(transPr1)) # update iterative parameters in an arrray
      # throw warning when any element in the iterative parameters can't be calculated

      ## 2. theta2 ---------------------------------------------------
      esAvg = esAvg1; eSize = eSize1; initPr = initPr1; transPr = transPr1
      Nmat = eSize%*%t(ones_nbT) # matrix of size
      Mmat = esAvg%*%t(ones_nbT) # matrix of mean
      Xmat = ones_nbS%*%t(X)     # matrix of sample X
      ## Estep --------
      Pmat = dnbinom(Xmat,size=Nmat,mu=Mmat) # matrix of emission Probability
      nbT = ncol(Pmat); nbS = nrow(Pmat) # forward
      Fpr = matrix(NA, nbS, nbT)
      Apr = rep(NA, nbT)
      Apr[1]  = sum(initPr*Pmat[,1])
      Fpr[,1]  = initPr*Pmat[,1]/Apr[1]
      for(tt in 2:nbT){
        Fpr.tmp = rep(0,nbS)
        for(l in 1:nbS){
          Fpr.tmp[l] = Pmat[l,tt]*sum(Fpr[,tt-1]*transPr[,l])
        }
        Apr[tt]  = sum(Fpr.tmp)
        Fpr[,tt] = Fpr.tmp/Apr[tt]
      }
      resFwd = list(Fpr=Fpr, logLik=sum(log(Apr)))
      Fpr = resFwd$Fpr
      logLik = resFwd$logLik
      nbT = ncol(Fpr); nbS = nrow(Fpr) # backward
      postPr = matrix(NA, nbS, nbT)
      Gpr = matrix(NA, nbS, nbT)
      postPr[,nbT] = Fpr[,nbT]
      for(tt in (nbT-1):1){
        Gpr[,tt+1]  = Fpr[,tt] %*% transPr
        postPr[,tt] = Fpr[,tt] *(transPr %*% (postPr[,tt+1]/Gpr[,tt+1]))#sapply(1:nbS, function(l) sum(trans[l,] * postPr[,tt+1] / Gpr[,tt+1]))#
      }
      resBcd  = list(Gpr = Gpr, postPr = postPr)
      postPr  = resBcd$postPr
      Gpr = resBcd$Gpr
      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
        esAvg = sort(eSize*(1-eProb)/eProb);  # new value of mean
        esVar = eSize*(1-eProb)/eProb^2;  # new value of variance
      }else if(algo=="GEM-NM"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          eParam = Nelder.Mead(objFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-BFGS"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = BFGS(objFunc, devFunc, par)$xmin
          # eParam = optim(par,objFunc,devFunc,method = "BFGS")$par
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-DFP"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = DFP(objFunc, devFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }
      Pi.tmp = matrix(NA, nbS, nbS)
      eParam = NULL
      for(k in 1:nbS){
        for(l in 1:nbS){
          Pi.tmp[k,l] = transPr[k,l]*sum( Fpr[k,1:(nbT-1)]*postPr[l,2:nbT]/Gpr[l,2:nbT] )
        }
        Mean.tmp = sum(postPr[k,] * X) / sum(postPr[k,])
        Var.tmp  = sum(postPr[k,] * (X - Mean.tmp)^2) / sum(postPr[k,])
        MV.tmp = c(Mean.tmp, Var.tmp)
        eParam = rbind(eParam, MV.tmp)
      }
      transPr = Pi.tmp/rowSums(Pi.tmp)
      initPr = postPr[,1]
      esAvg2 = esAvg; eSize2 = eSize; initPr2 = initPr; transPr2 = transPr
      newParam2 = c(esAvg2, eSize2, initPr2, c(transPr2)) # update iterative parameters in an arrray
      # throw warning when any element in the iterative parameters can't be calculated

      ## 3. r ---------------------------------------------------
      r = newParam1 - prevParam

      ## 4. v ---------------------------------------------------
      v = newParam2 - 2*newParam1 + prevParam

      ## 5. alpha  ---------------------------------------------------
      if(step==1){
        alpha= sum(r*v) / sum(v*v)
      }else if (step==2){
        alpha= sum(r*r) / sum(r*v)
      }else if (step==3){
        alpha = - sqrt(sum(r * r))/ sqrt(sum(v * v))
      }else{
        warning("Please check the 'step' input!")
      }

      ## 6. theta prim ---------------------------------------------------
      thetaprim = prevParam -2*alpha*r + alpha^2*v
      esAvgprim = thetaprim[1:nbS]
      eSizeprim = thetaprim[(nbS+1):(2*nbS)]
      initPrprim = thetaprim[(2*nbS+1):(3*nbS)]
      transPrprim = thetaprim[(3*nbS+1):(4*nbS)]
      newParamprim = c(esAvgprim, eSizeprim, initPrprim, transPrprim)

      ## 7. theta0 ---------------------------------------------------
      esAvg = esAvgprim; eSize = eSizeprim; initPr = initPrprim; transPr = matrix(c(transPrprim),nbS,nbS)
      Nmat = eSize%*%t(ones_nbT) # matrix of size
      Mmat = esAvg%*%t(ones_nbT) # matrix of mean
      Xmat = ones_nbS%*%t(X)     # matrix of sample X
      ## Estep --------
      Pmat = dnbinom(Xmat,size=Nmat,mu=Mmat) # matrix of emission Probability
      nbT = ncol(Pmat); nbS = nrow(Pmat) # forward
      Fpr = matrix(NA, nbS, nbT)
      Apr = rep(NA, nbT)
      Apr[1]  = sum(initPr*Pmat[,1])
      Fpr[,1]  = initPr*Pmat[,1]/Apr[1]
      for(tt in 2:nbT){
        Fpr.tmp = rep(0,nbS)
        for(l in 1:nbS){
          Fpr.tmp[l] = Pmat[l,tt]*sum(Fpr[,tt-1]*transPr[,l])
        }
        Apr[tt]  = sum(Fpr.tmp)
        Fpr[,tt] = Fpr.tmp/Apr[tt]
      }
      resFwd = list(Fpr=Fpr, logLik=sum(log(Apr)))
      Fpr = resFwd$Fpr
      logLik = resFwd$logLik
      nbT = ncol(Fpr); nbS = nrow(Fpr) # backward
      postPr = matrix(NA, nbS, nbT)
      Gpr = matrix(NA, nbS, nbT)
      postPr[,nbT] = Fpr[,nbT]
      for(tt in (nbT-1):1){
        Gpr[,tt+1]  = Fpr[,tt] %*% transPr
        postPr[,tt] = Fpr[,tt] *(transPr %*% (postPr[,tt+1]/Gpr[,tt+1]))#sapply(1:nbS, function(l) sum(trans[l,] * postPr[,tt+1] / Gpr[,tt+1]))#
      }
      resBcd  = list(Gpr = Gpr, postPr = postPr)
      postPr  = resBcd$postPr
      Gpr = resBcd$Gpr
      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
        esAvg = sort(eSize*(1-eProb)/eProb);  # new value of mean
        esVar = eSize*(1-eProb)/eProb^2;  # new value of variance
      }else if(algo=="GEM-NM"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          eParam = Nelder.Mead(objFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-BFGS"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = BFGS(objFunc, devFunc, par)$xmin
          # eParam = optim(par,objFunc,devFunc,method = "BFGS")$par
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }else if(algo=="GEM-DFP"){
        ## Mstep ————
        for(k in 1:nbS){
          par = c(eSize[k], esAvg[k])
          objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
          devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
          eParam = DFP(objFunc, devFunc, par)$xmin
          eSize[k] = eParam[1]; esAvg[k] = eParam[2]
        }
      }
      Pi.tmp = matrix(NA, nbS, nbS)
      eParam = NULL
      for(k in 1:nbS){
        for(l in 1:nbS){
          Pi.tmp[k,l] = transPr[k,l]*sum( Fpr[k,1:(nbT-1)]*postPr[l,2:nbT]/Gpr[l,2:nbT] )
        }
        Mean.tmp = sum(postPr[k,] * X) / sum(postPr[k,])
        Var.tmp  = sum(postPr[k,] * (X - Mean.tmp)^2) / sum(postPr[k,])
        MV.tmp = c(Mean.tmp, Var.tmp)
        eParam = rbind(eParam, MV.tmp)
      }
      transPr = Pi.tmp/rowSums(Pi.tmp)
      initPr = postPr[,1]
      newParam = c(esAvg, eSize, initPr, c(transPr)) # update iterative parameters in an arrray
      # throw warning when any element in the iterative parameters can't be calculated

      ## 8. convergence ---------------------------------------------------
      if(any(is.na(esAvg)) || any(is.na(esAvg))){
        warning("HMM EM - NA is appeared!")
        esAvg=prevParam[1:nbS]; eSize=prevParam[(nbS+1):(2*nbS)];
        break
      }
      D = sum((newParam - prevParam)^2) # when precision reach the desired level, stop iteration, otherwise reach the maximum limit
      if(D < EPSILON){break}
    }
    cat("iterStop = ", iter, "\n")
    # return with a list
    parameter=list(); parameter[["eSize"]]=eSize; parameter[["eProb"]]=eProb; parameter[["esAvg"]]=esAvg;
    parameter[["esVar"]]=esVar; parameter[["initPr"]]=as.numeric(initPr); parameter[["transPr"]]=as.numeric(transPr)
    structure(list(parameter=parameter, emisPr=Pmat, postPr=postPr), class="HMMNB")

  }else if(accelerate=="DAAREM"){
    if(algo=="EM-FL"){
      para.init <- c(esAvg, eSize, initPr, c(transPr))
      res<-HMMNB.daarem(par=para.init, fixptfn = HMMNB.EMFL.Update, objfn = HMMNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))
      return(res)

    }else if(algo=="GEM-NM"){
      para.init <- c(esAvg, eSize, initPr, c(transPr))
      res<-HMMNB.daarem(par=para.init, fixptfn = HMMNB.NM.Update, objfn = HMMNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))
      return(res)
    }else if(algo=="GEM-BFGS"){
      para.init <- c(esAvg, eSize, initPr, c(transPr))
      res<-HMMNB.daarem(par=para.init, fixptfn = HMMNB.BFGS.Update, objfn = HMMNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))
      return(res)
    }else if(algo=="GEM-DFP"){
      para.init <- c(esAvg, eSize, initPr, c(transPr))
      res<-HMMNB.daarem(par=para.init, fixptfn = HMMNB.DFP.Update, objfn = HMMNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))
      return(res)
    }

  }


}

