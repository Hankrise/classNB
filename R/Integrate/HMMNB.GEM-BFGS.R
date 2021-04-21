HMMNB.BFGS.Update <- function(para.hat, X, nbS) {
  nbT=length(X)
  esAvg<-para.hat[1:nbS]
  eSize<-para.hat[(nbS+1):(2*nbS)]
  initPr<-para.hat[(2*nbS+1):(3*nbS)]
  transPr<-matrix(para.hat[(3*nbS+1):((nbS)*(nbS+2))],nbS,nbS,byrow=T)
  esVar = esAvg^2/eSize+esAvg # absence cause error
  eProb = esAvg/esVar
  ones_nbS=as.matrix(rep(1,nbS))
  ones_nbT=as.matrix(rep(1,nbT))
  ones_nbTS=ones_nbT%*%t(ones_nbS)
  ones_nbST=ones_nbS%*%t(ones_nbT)
  ones_nbSS=ones_nbS%*%t(ones_nbS)
  Nmat=eSize%*%t(ones_nbT) # size
  Mmat=esAvg%*%t(ones_nbT) # avg
  Xmat=ones_nbS%*%t(X)     # X
  ## Estep --------
  Pmat=dnbinom(Xmat,size = Nmat, mu=Mmat) # emisPr
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
  beta = 1 - 1/(1-eProb) - 1/log(eProb) # value of beta
  Bmat = beta%*%t(ones_nbT) # matrix of beta
  Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #matrix of delta
  ## Mstep ————
  for(k in 1:nbS){
    par = c(eSize[k], esAvg[k])
    objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
    devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
    eParam = BFGS(objFunc, devFunc, par)$xmin
    eSize[k] = eParam[1]; esAvg[k] = eParam[2]
  }
  sumWeight = rowSums(postPr); eWeight = sumWeight/sum(sumWeight)
  # esAvg = eSize * (1-eProb)/eProb;
  # esVar = eSize * (1-eProb)/eProb^2;
  eProb = eSize/(eSize+esAvg)
  esVar = esAvg^2/eSize+esAvg
  ordAvg = order(esAvg)
  esAvg = esAvg[ordAvg]
  eSize = eSize[ordAvg]
  eProb = eProb[ordAvg]
  esVar = esVar
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
  newParam = c(esAvg, eSize, initPr, c(transPr))
  return(newParam)
}
