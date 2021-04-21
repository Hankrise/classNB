#' Parameter initialisation 
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of mean (optional).
#' @param var a vector of variance (optional).
#' @return A list of 5 objets.
#' \describe{
#' \item{\code{eWeight}}{a vector of the initialised parameter of weight.}
#' \item{\code{esAvg}}{a vector of the initialised parameter of mean.}
#' \item{\code{esVar}}{a vector of the initialised parameter of variance.}
#' \item{\code{eSize}}{a vector of the initialised parameter of size.}
#' \item{\code{eProb}}{a vector of the initialised parameter of probability.}
#' }
mixNB.init <- function(X, nbS,avg=NULL, var=NULL){
  eWeight = rep(1, nbS)/nbS # the weight in NB mixture model
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

  return(list(eWeight=eWeight,esAvg=esAvg,esVar=esVar,eSize=eSize,eProb=eProb))
}
#' Perform inference of mixture of negative binomial models
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of initial value of mean (optional).
#' @param var a vector of initial value of variance (optional).
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_CF algorithm (optional).
#' @param EPSILON a value for the threshold used for the stopping criteria for the mixNB_CF algorithm (optional).
#' @param algo a string of specifying the approach of M-step in EM algorithm, available choices are "EM-FL", "GEM-NM", "GEM-BFGS", "GEM-DFP".
#' @param accelerate a string of specifying the acceleration method of EM algorithm, available choices are NULL(no accleration), "SQEM", "DAAREM".
#' @return A list of 3 objets.
#' \describe{
#' \item{\code{parameter}}{a list with the inferred parameters: Size, Probability, mean, variance, weight.}
#' \item{\code{emisPr}}{a matrix with emission probability associated to each states in column and each position in row.}
#' \item{\code{postPr}}{a matrix with posterior probability associated to each series in column and each position in row.}
#' \item{\code{iter}}{a number of iterations while algorithm stop.}
#' }
#' @export
#' @examples
#' data(toyexample)
#' # Perform inference of mixture of negative binomial models
#' resmixNB <- mixNB(X = toydata, nbS = 3,algo="EM-FL")
#' @seealso \code{\link{mixNB_NS}}, \code{\link{mixNB_NS}}
mixNB <- function(X, nbS,avg=NULL, var=NULL, iterMax=5000, EPSILON=1e-7,algo="EM-FL",accelerate=NULL){
  convergence<-0
  init=mixNB.init(X, nbS,avg, var)
  eWeight=init$eWeight
  nbT=length(X)
  esAvg=init$esAvg
  esVar=init$esVar
  eSize=init$eSize
  eProb=init$eProb
  ones_nbS=as.matrix(rep(1,nbS))
  ones_nbT=as.matrix(rep(1,nbT))
  ones_nbTS=ones_nbT%*%t(ones_nbS)
  ones_nbST=ones_nbS%*%t(ones_nbT)
  ones_nbSS=ones_nbS%*%t(ones_nbS)
  if(is.null(accelerate)){
    for(iter in 1:iterMax){
      prevParam = c(esAvg, eSize, eWeight)
      Nmat=eSize%*%t(ones_nbT) # size
      Mmat=esAvg%*%t(ones_nbT) # avg
      Xmat=ones_nbS%*%t(X)     # X
      ## Estep --------
      Pmat=dnbinom(Xmat,size = Nmat, mu=Mmat) # emisPr
      rawPost = Pmat*eWeight; postPr = rawPost/(ones_nbSS%*%rawPost)
      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
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
      sumWeight=diag(postPr%*%ones_nbTS);eWeight=sumWeight/sum(sumWeight)
      if(algo=="EM-FL"){
        esAvg = eSize * (1-eProb)/eProb;
        esVar = eSize * (1-eProb)/eProb^2;
      }
      ordAvg = order(esAvg)
      esAvg = esAvg[ordAvg]
      eSize = eSize[ordAvg]
      eWeight= eWeight[ordAvg]
      if(algo=="EM-FL"){
        eProb = eProb[ordAvg]
        esVar = esVar[ordAvg]
      }else{
        eProb = eSize/(eSize+esAvg)
        esVar = esAvg^2/eSize+esAvg
      }
      newParam = c(esAvg, eSize, eWeight)
      if(any(is.na(esAvg)) || any(is.na(eSize))){
        # warning("mixNB - NA is appeared!")
        convergence<-0
        esAvg=prevParam[1:nbS]; eSize=prevParam[(nbS+1):(2*nbS)];
        break
      }
      D = sum((newParam - prevParam)^2)
      if(D < EPSILON){  convergence<-1;break}
    }
    cat("iterStop = ", iter, "\n")
    parameter=list(); parameter[["eSize"]]=eSize; parameter[["eProb"]]=eProb; parameter[["esAvg"]]=esAvg;
    parameter[["esVar"]]=esVar; parameter[["eWeight"]]=eWeight;
    structure(list(parameter=parameter, emisPr=Pmat, postPr=postPr,iter=iter,convergence=convergence), class="mixNB")


  }else if(accelerate=="SQEM"){
    step=1
    for(iter in 1:iterMax){
      prevParam = c(esAvg, eSize, eWeight) # define iterative parameters in an arrray

      ## 1. theta1 ---------------------------------------------------
      Nmat=eSize%*%t(ones_nbT) # matrix of size
      Mmat=esAvg%*%t(ones_nbT) # matrix of mean
      Xmat=ones_nbS%*%t(X)     # matrix of sample X
      # Estep --------
      Pmat=dnbinom(Xmat,size = Nmat, mu=Mmat) # matrix of emission Probability
      rawPost = Pmat*eWeight; postPr = rawPost/(ones_nbSS%*%rawPost)# matrix of post Probability
      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
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
      sumWeight=diag(postPr%*%ones_nbTS);eWeight=sumWeight/sum(sumWeight) # new value of weight
      if(algo=="EM-FL"){
        esAvg = eSize * (1-eProb)/eProb;
        esVar = eSize * (1-eProb)/eProb^2;
      }

      ordAvg = order(esAvg)
      esAvg = esAvg[ordAvg]
      eSize = eSize[ordAvg]
      eWeight= eWeight[ordAvg]
      if(algo=="EM-FL"){
        eProb = eProb[ordAvg]
        esVar = esVar[ordAvg]
      }else{
        eProb = eSize/(eSize+esAvg)
        esVar = esAvg^2/eSize+esAvg
      }
      esAvg1 = esAvg; eSize1 = eSize; eWeight1 = eWeight
      newParam1 = c(esAvg1, eSize1, eWeight1) # update iterative parameters in an arrray
      # throw warning when any element in the iterative parameters can't be calculated

      ## 2. theta2 ---------------------------------------------------
      esAvg = esAvg1; eSize = eSize1; eWeight = eWeight1;
      Nmat=eSize%*%t(ones_nbT) # matrix of size
      Mmat=esAvg%*%t(ones_nbT) # matrix of mean
      Xmat=ones_nbS%*%t(X)     # matrix of sample X
      # Estep --------
      Pmat=dnbinom(Xmat,size = Nmat, mu=Mmat) # matrix of emission Probability
      rawPost = Pmat*eWeight; postPr = rawPost/(ones_nbSS%*%rawPost)# matrix of post Probability
      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
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
      sumWeight=diag(postPr%*%ones_nbTS);eWeight=sumWeight/sum(sumWeight) # new value of weight
      if(algo=="EM-FL"){
        esAvg = eSize * (1-eProb)/eProb;
        esVar = eSize * (1-eProb)/eProb^2;
      }

      ordAvg = order(esAvg)
      esAvg = esAvg[ordAvg]
      eSize = eSize[ordAvg]
      eWeight= eWeight[ordAvg]
      if(algo=="EM-FL"){
        eProb = eProb[ordAvg]
        esVar = esVar[ordAvg]
      }else{
        eProb = eSize/(eSize+esAvg)
        esVar = esAvg^2/eSize+esAvg
      }
      esAvg2 = esAvg; eSize2 = eSize; eWeight2 = eWeight;
      newParam2 = c(esAvg2, eSize2, eWeight2) # update iterative parameters in an arrray
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
      eWeightprim = thetaprim[(2*nbS+1):(3*nbS)]
      newParamprim = c(esAvgprim, eSizeprim, eWeightprim)

      ## 7. theta0 ---------------------------------------------------
      esAvg = esAvgprim; eSize = eSizeprim; eWeight = eWeightprim;
      Nmat=eSize%*%t(ones_nbT) # matrix of size
      Mmat=esAvg%*%t(ones_nbT) # matrix of mean
      Xmat=ones_nbS%*%t(X)     # matrix of sample X
      # Estep --------
      Pmat=dnbinom(Xmat,size = Nmat, mu=Mmat) # matrix of emission Probability
      rawPost = Pmat*eWeight; postPr = rawPost/(ones_nbSS%*%rawPost)# matrix of post Probability

      if(algo=="EM-FL"){
        beta = 1 - 1/(1-eProb) - 1/log(eProb)
        Bmat=beta%*%t(ones_nbT) # beta
        Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
        ## Mstep ---------
        eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
        eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
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
      sumWeight=diag(postPr%*%ones_nbTS);eWeight=sumWeight/sum(sumWeight) # new value of weight
      if(algo=="EM-FL"){
        esAvg = eSize * (1-eProb)/eProb;
        esVar = eSize * (1-eProb)/eProb^2;
      }

      ordAvg = order(esAvg)
      esAvg = esAvg[ordAvg]
      eSize = eSize[ordAvg]
      eWeight= eWeight[ordAvg]
      if(algo=="EM-FL"){
        eProb = eProb[ordAvg]
        esVar = esVar[ordAvg]
      }else{
        eProb = eSize/(eSize+esAvg)
        esVar = esAvg^2/eSize+esAvg
      }
      newParam = c(esAvg, eSize, eWeight) # update iterative parameters in an arrray
      # throw warning when any element in the iterative parameters can't be calculated

      ## 8. convergence ---------------------------------------------------
      if(any(is.na(esAvg)) || any(is.na(esAvg))){
        # warning("mixNB - NA is appeared!")
        convergence<-0
        esAvg=prevParam[1:nbS]; eSize=prevParam[(nbS+1):(2*nbS)];
        break
      }
      D = sum((newParam - prevParam)^2) # when precision reach the desired level, stop iteration, otherwise reach the maximum limit
      if(D < EPSILON){  convergence<-1;break}
    }
    cat("iterStop = ", iter, "\n")
    # return with a list
    parameter=list(); parameter[["eSize"]]=eSize; parameter[["eProb"]]=eProb; parameter[["esAvg"]]=esAvg;
    parameter[["esVar"]]=esVar; parameter[["eWeight"]]=eWeight;
    structure(list(parameter=parameter, emisPr=Pmat, postPr=postPr,iter=iter,convergence=convergence), class="mixNB")

  }else if(accelerate=="DAAREM"){
    if(algo=="EM-FL"){
      para.init <- c(esAvg,eSize,eWeight)
      res<-mixNB.daarem(par=para.init, fixptfn = mixNB.EMFL.Update, objfn = mixNBLogLik,
                           X=X, nbS=nbS, control=list(maxiter=iterMax))
      
    }else if(algo=="GEM-NM"){
      para.init <- c(esAvg,eSize,eWeight)
      res<-mixNB.daarem(par=para.init, fixptfn = mixNB.NM.Update, objfn = mixNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))

    }else if(algo=="GEM-BFGS"){
      para.init <- c(esAvg,eSize,eWeight)
      res<-mixNB.daarem(par=para.init, fixptfn = mixNB.BFGS.Update, objfn = mixNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))

    }else if(algo=="GEM-DFP"){
      para.init <- c(esAvg,eSize,eWeight)
      res<-mixNB.daarem(par=para.init, fixptfn = mixNB.DFP.Update, objfn = mixNBLogLik,
                        X=X, nbS=nbS, control=list(maxiter=iterMax))

    }
    
    esAvg=res$par[1:nbS]
    eSize=res$par[(1:nbS)+nbS]
    eWeight=res$par[(1:nbS)+nbS*2]
    eProb=eSize/(eSize+esAvg)
    esVar=esAvg^2/eSize+esAvg
    parameter=list(); parameter[["eSize"]]=eSize; parameter[["eProb"]]=eProb; parameter[["esAvg"]]=esAvg;
    parameter[["esVar"]]=esVar; parameter[["eWeight"]]=eWeight;
    ones_nbS=as.matrix(rep(1,nbS))
    ones_nbT=as.matrix(rep(1,nbT))
    ones_nbSS=ones_nbS%*%t(ones_nbS)
    Nmat=eSize%*%t(ones_nbT) # matrix of size
    Mmat=esAvg%*%t(ones_nbT) # matrix of mean
    Xmat=ones_nbS%*%t(X)     # matrix of sample X
    Pmat=dnbinom(Xmat,size = Nmat, mu=Mmat)
    rawPost = Pmat*eWeight; postPr = rawPost/(ones_nbSS%*%rawPost)
    iter=length(res[[6]])
    convergence=res$convergence
    structure(list(parameter=parameter, emisPr=Pmat, postPr=postPr,iter=iter,convergence=convergence), class="mixNB")
    }


  }

