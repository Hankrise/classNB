#' Initialisation HMM 
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param meth.init a string specifying the method of initialization.
#' @param avg a vector of initial value of mean (optional).
#' @param var a vector of initial value of variance (optional).
#' @return A list of 3 objets.
#' \describe{
#' \item{\code{init}}{a vector of initial value of initial probability of HMM.}
#' \item{\code{trans}}{a matrix of initial value of transition probability matrix of HMM.}
#' \item{\code{eParam}}{a matrix of the inferred parameter: average and variance bind into a matrix.}
#' }
initHMM <- function(X, nbS, meth.init="kmeans",avg=NULL,var=NULL)
{
  init  = rep(1/nbS, nbS)
  trans = 0.5*diag(nbS) + array(0.5/(nbS),c(nbS,nbS))
  
  if(is.null(c(avg,var))){
    # if initial value of mean & variance are not given, use kmeans algorithm to generate
    esAvg=1;esVar=0; # ensure mean < variance
    while(any(esAvg>esVar)){
      clust <- kmeans(x = as.vector(t(X)), centers = nbS) # apply kmeans
      esAvg = as.vector(clust$centers) # calculate initial mean 
      ordAvg = order(esAvg) # the order of initial mean
      esAvg = esAvg[ordAvg] # rearrange mean in increasing order
      esVar = clust$withinss / (clust$size - 1) # calculate initial variance
      esVar = esVar[ordAvg] # rearrange variance in increasing order of mean
    }
  } else {
    # if initial value is given, use them.
    esAvg = avg
    esVar = var
  }
  
  return(list(init=init, trans=trans, eParam=cbind(esAvg, esVar)))
}

#' Calculate emission probability
#' @param X a vector of observations.
#' @param eParam a matrix of the inferred parameter: average and variance bind into a matrix.
#' @return A matrix of size K*T.
#' \describe{
#' \item{\code{emisPr}}{an emission distribution matrix.}
#' }
EmisPr <- function(X, eParam)
{ #eParam : a matrix of size K*2 storing the mean and variance
  nbS   = nrow(eParam); nbT=length(X)
  emisPr = matrix(NA, nbS, nbT)
  for(t in 1:nbT){
    emisPr[,t]  = sapply(1:nbS, function(k) dnorm(X[t], eParam[k,1], sqrt(eParam[k,2])))  
  }
   emisPr
} 

#' Calculate forward probability of HMM
#' @param init a vector of initial probability of HMM.
#' @param trans a matrix of transition probability matrix of HMM.
#' @param emisPr a matrix of emission probability.
#' @return A list of 2 objets.
#' \describe{
#' \item{\code{Fpr}}{a matrix of forward probability of HMM in log scale.}
#' \item{\code{logLik}}{a numeric value of the log liklihood value.}
#' }
Forward <- function(init, trans, emisPr)
{
  nbT = ncol(emisPr); nbS = nrow(emisPr)
  Fpr    = matrix(NA, nbS, nbT)
  Apr    = rep(NA, nbT)
  
   Apr[1]  = sum(init*emisPr[,1])
  Fpr[,1]  = init*emisPr[,1]/Apr[1]
  for(tt in 2:nbT){
    Fpr.tmp = rep(0,nbS)
    for(l in 1:nbS){
      Fpr.tmp[l] = emisPr[l,tt]*sum(Fpr[,tt-1]*trans[,l])
    }
    Apr[tt]  = sum(Fpr.tmp)
    Fpr[,tt] = Fpr.tmp/Apr[tt]
  }  
  # print(Apr)
  return(list(Fpr=Fpr, logLik=sum(log(Apr))))
} 

#' Calculate backward probability of HMM
#' @param Fpr a matrix of forward probability of HMM.
#' @param trans a matrix of transition probability matrix of HMM.
#' @return A list of 2 objets.
#' \describe{
#' \item{\code{Gpr}}{a matrix of backward probability of HMM.}
#' \item{\code{postPr}}{a matrix of posterior probability.}
#' }
Backward <- function(Fpr, trans)
{
  nbT = ncol(Fpr); nbS = nrow(Fpr)
  postPr    = matrix(NA, nbS, nbT)
     Gpr    = matrix(NA, nbS, nbT)
  
  postPr[,nbT] = Fpr[,nbT]
  for(tt in (nbT-1):1){
      Gpr[,tt+1]  = Fpr[,tt] %*% trans
      postPr[,tt] = Fpr[,tt] *(trans %*% (postPr[,tt+1]/Gpr[,tt+1]))#sapply(1:nbS, function(l) sum(trans[l,] * postPr[,tt+1] / Gpr[,tt+1]))#  
  }
  # print(t(round(postPr,2)))
  return(list(Gpr = Gpr, postPr = postPr))
} 

#' E-step of expectation maximization
#' @param init a vector of initial probability of HMM.
#' @param trans a matrix of transition probability matrix of HMM.
#' @param emisPr a matrix of emission probability.
#' @param eParam a matrix of the inferred parameter: average and variance bind into a matrix.
#' @return A list of 5 objets.
#' \describe{
#' \item{\code{Fpr}}{a matrix of forward probability of HMM in log scale.}
#' \item{\code{Gpr}}{a matrix of backward probability of HMM.}
#' \item{\code{postPr}}{a matrix of posterior probability.}
#' \item{\code{logLik}}{a numeric value of the log liklihood value.}
#' \item{\code{eParam}}{a matrix of the inferred parameter: average and variance bind into a matrix.}
#' }
Estep <- function(init, trans, emisPr,eParam)
{
  nbT = ncol(emisPr); nbS = nrow(emisPr)
  
  resFwd  <- Forward(init, trans, emisPr)
       Fpr   = resFwd$Fpr
    logLik   = resFwd$logLik
  resBcd  <- Backward(Fpr, trans)
    postPr   = resBcd$postPr
       Gpr   = resBcd$Gpr
  
  return(list(Fpr=Fpr, Gpr=Gpr, postPr=postPr, logLik=logLik,eParam=eParam))
} 

#' M-step of expectation maximization with numeric solution
#' @param Fpr a matrix of forward probability of HMM in log scale.
#' @param Gpr a matrix of backward probability of HMM.
#' @param postPr a matrix of posterior probability.
#' @param trans a matrix of transition probability matrix of HMM.
#' @param X a vector of observations.
#' @param eParam a matrix of the inferred parameter: average and variance bind into a matrix.
#' @return A list of 5 objets.
#' \describe{
#' \item{\code{init}}{a vector of the undated initial probability of HMM.}
#' \item{\code{trans}}{a matrix of the undated transition probability matrix of HMM.}
#' \item{\code{eParam}}{a matrix of the undated inferred parameter: average and variance bind into a matrix.}
#' }
Mstep <- function(Fpr, Gpr, postPr, trans, X,eParam)
{
  nbT = ncol(Fpr); nbS = nrow(Fpr)
  Pi.tmp = matrix(NA, nbS, nbS)
  # eParam = NULL
  for(k in 1:nbS){
    for(l in 1:nbS){
      Pi.tmp[k,l] = trans[k,l]*sum( Fpr[k,1:(nbT-1)]*postPr[l,2:nbT]/Gpr[l,2:nbT] )
    } 
    # Mean.tmp = sum(postPr[k,] * X) / sum(postPr[k,])
    # Var.tmp  = sum(postPr[k,] * (X - Mean.tmp)^2) / sum(postPr[k,])
    # MV.tmp = c(Mean.tmp, Var.tmp)
    # eParam = rbind(eParam, MV.tmp)
  }
  for(k in 1:nbS){
    par = c(eParam[k,])
    objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2]) * postPr[k,])}
    # devFunc <- function(par){-as.matrix(c(sum(postPr[k,]*(digamma(X+par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k,]*(par[1]/par[2] - X/(1-par[2])))))}
    # eParam = BFGS(objFunc, devFunc, par)$xmin
    # newpar = optim(par,objFunc,devFunc,method = "BFGS")$par
    newpar = optim(par,objFunc,method = "BFGS")$par
    # eSize[k] = eParam[1]; esAvg[k] = eParam[2]
    eParam[k,]=newpar
  }
  trans = Pi.tmp/rowSums(Pi.tmp)
  # ordMean = order(eParam[,1])
  # init = postPr[,1][ordMean]
  # eParam = eParam[ordMean,]
  # trans = trans[ordMean, ordMean]
  init = postPr[,1]

  return(list(init=init, trans=trans, eParam=eParam))
} 

#' Expectation maximization with numeric approach in M-step
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_CF algorithm.
#' @param threshold a value for the threshold used for the stopping criteria for the mixNB_CF algorithm.
#' @return A list of 5 objets.
#' \describe{
#' \item{\code{initPr}}{a vector of the undated initial probability of HMM.}
#' \item{\code{transPr}}{a matrix of the undated transition probability matrix of HMM.}
#' \item{\code{eParam}}{a matrix of the undated inferred parameter: average and variance bind into a matrix.}
#' \item{\code{postPr}}{a matrix of posterior probability.}
#' \item{\code{logLik}}{a numeric value of the log liklihood value.}
#' \item{\code{emisPr}}{an emission distribution matrix.}
#' \item{\code{iter}}{a number of iterations while algorithm stop.}
#' }
EM <- function(X, nbS, iterMax, threshold)
{
  resInitHMM <- initHMM(X,nbS)
    init    = resInitHMM$init
    trans   = resInitHMM$trans
    eParam  = resInitHMM$eParam
    Param   = c(init, as.vector(trans), as.vector(eParam))

  for(iter in 1:iterMax){
    emisPr <- EmisPr(X, eParam)
    
    resEstep <- Estep(init, trans, emisPr,eParam)
         Fpr  = resEstep$Fpr
         Gpr  = resEstep$Gpr
      postPr  = resEstep$postPr
      logLik  = resEstep$logLik
    
    resMstep <- Mstep(Fpr, Gpr, postPr, trans, X,eParam)
      init.new   = resMstep$init
      trans.new  = resMstep$trans
      eParam.new = resMstep$eParam
    
    Param.new = c(init.new, as.vector(trans.new), as.vector(eParam.new))  
    D = sum((Param.new - Param)^2)
    Param  = Param.new
    init   = init.new
    trans  = trans.new
    eParam = eParam.new
    
    eParam = eParam[order(eParam[,1]),]

    if(D<threshold){break}
  }
  return(list(initPr=init, transPr=trans, eParam=eParam, postPr=postPr, logLik=logLik, emisPr=emisPr, iter=iter))
} 


#' Expectation maximization with closed-form approach in M-step
#' @param initPr a vector of the initial probability of HMM.
#' @param transPr a matrix of the transition probability matrix of HMM.
#' @param emisPr amatrix of emission probability.
#' @return A numeric vector
#' \describe{
#' \item{\code{path}}{a vector of the inferred states of each observations.}
#' }
Viterbi <- function(initPr,transPr,emisPr)
{
  nbT = ncol(emisPr); nbS = nrow(emisPr)
  logF = matrix(0, ncol=nbT,nrow=nbS)
  logF[,1] = log(initPr) + log(emisPr[,1])
  
  for (tt in 2:nbT) {
    for (k in 1:nbS) {
      logF[k,tt] = log(emisPr[k,tt]) +  max(logF[,tt-1] + log(transPr[,k]))
    }
  }
  path = rep(NA, nbT)
  path[nbT] = which.max(logF[,nbT])
  
  for (tt in (nbT-1):1) {
    path[tt] = which.max(logF[,tt] + log(transPr[,path[tt+1]]))
  }
  return(path)
}

#' Perform inference of HMM negative binomial model using numeric approach
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param viterbi a logical variable indicating whether to use Maximum A Posteriori method (FALSE) or Viterbi algorithm (TRUE, by default) for recovering the most likely path.
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_CF algorithm (optional).
#' @param threshold a value for the threshold used for the stopping criteria for the mixNB_CF algorithm (optional).
#' @return A list of 8 objets.
#' \describe{
#' \item{\code{initPr}}{a vector of the inferred parameter: probability of initial distribution of HMM.}
#' \item{\code{transPr}}{a matrix of the inferred parameter: transition probability matrix.}
#' \item{\code{eParam}}{a matrix of the inferred parameter: average and variance bind into a matrix.}
#' \item{\code{postPr}}{a matrix of posterior probability.}
#' \item{\code{logLik}}{a numeric value of the log liklihood value.}
#' \item{\code{emisPr}}{a matrix of emission probability.}
#' \item{\code{iterStop}}{a number of iterations while algorithm stop.}
#' \item{\code{path}}{a vector of the inferred states of each observations.}
#' }
#' @export
#' @examples
#' data(toyexample)
#' # Perform inference of mixture of negative binomial models using numeric approach
#' resHMM_NS <- HMM_NS(X = toydata, nbS = 3)
#' @seealso \code{\link{HMM_CF}}, \code{\link{HMM_CF}}
HMM_NB_NS <- function(X, nbS, viterbi=TRUE,iterMax=500, threshold=1e-7){
  resEM <- EM(X, nbS, iterMax, threshold)
  postPr = resEM$postPr
  initPr = resEM$initPr
  transPr = resEM$transPr
  eParam  = resEM$eParam
  emisPr  =resEM$emisPr
  iterStop = resEM$iter
  logLik  = resEM$logLik
  if(viterbi){
    path <- Viterbi(initPr,transPr,emisPr)
  }else{
    path = apply(postPr,2,which.max)
  }
  return(list(initPr=initPr, transPr=transPr, eParam=eParam, postPr=postPr, logLik=logLik, emisPr=emisPr,iterStop=iterStop, path=path))
}

