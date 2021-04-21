# EM-FL iter on eProb->esAvg, eSize, eWeight
mixNB.EMFL.Update <- function(para.hat, X, nbS) {
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
  ## Mstep ---------
  eProb=beta*diag(postPr%*%t(Dmat))/diag(postPr%*%t(Xmat-((ones_nbST-Bmat)*Dmat)))
  eSize=-diag(postPr%*%t(Dmat))/diag(postPr%*%ones_nbTS)/log(eProb)
  sumWeight=diag(postPr%*%ones_nbTS);eWeight=sumWeight/sum(sumWeight)
  esAvg = eSize * (1-eProb)/eProb;
  esVar = eSize * (1-eProb)/eProb^2;
  ordAvg = order(esAvg)
  esAvg = esAvg[ordAvg]
  eSize = eSize[ordAvg]
  eWeight= eWeight[ordAvg]
  eProb = eProb[ordAvg]
  esVar = esVar[ordAvg]
  newParam = c(esAvg, eSize, eWeight)
  return(newParam)
}
