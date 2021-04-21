# GEM-DFP iter on esAvg, eSize, eWeight
mixNB.DFP.Update <- function(para.hat, X, nbS) {
  nbT=length(X)
  esAvg<-para.hat[1:nbS]
  eSize<-para.hat[(nbS+1):(2*nbS)]
  eWeight<-para.hat[(2*nbS+1):(3*nbS)]
  esVar = esAvg^2/eSize+esAvg
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
  rawPost = Pmat*eWeight; postPr = rawPost/(ones_nbSS%*%rawPost)
  beta = 1 - 1/(1-eProb) - 1/log(eProb)
  Bmat=beta%*%t(ones_nbT) # beta
  Dmat=Nmat*(digamma(Nmat+Xmat)-digamma(Nmat)) #delta
  ## Mstep ————
  for(k in 1:nbS){
    par = c(eSize[k], esAvg[k])
    objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
    devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
    eParam = DFP(objFunc, devFunc, par)$xmin
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
  eWeight= eWeight[ordAvg]
  eProb = eProb[ordAvg]
  esVar = esVar[ordAvg]
  newParam = c(esAvg, eSize, eWeight)
  return(newParam)
}
