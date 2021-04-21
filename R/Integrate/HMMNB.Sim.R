## -------------------------- ##
## Generate Markov chain
## leng   :: length of Markov chain
## trans  :: transition probability matrix
## init   :: initial distribution
## stationary :: TRUE/FALSE
## -------------------------- ##

MC <- function(leng, trans, init=NULL)
{ ## Set line state
  if(is.null(init)==TRUE){
    val.propre  <-  round(eigen(t(trans))$values,3)
    pos         <-  which(val.propre == 1.000)
    init.tmp    <-  eigen(t(trans))$vectors[,pos]
    init        <-  init.tmp / sum(init.tmp)
    init        <-  as.numeric(init)
  }else{
    init=init
  }

  state <- numeric(leng)
  state[1] <- which(rmultinom(1, size = 1, prob = init)==1)

  for ( tt in 2:leng){
    proba =  trans[state[tt-1], ]
    state[tt] <- which(rmultinom(1, size = 1, prob = proba )==1)
  }

  return(state)
}




HMMNB.Sim <- function(leng, trans, eMean, eVar, init=NULL,seed=2020)
{
  mc = MC(leng,trans,init)
  set.seed(seed)
  nbS=ncol(trans)
  eSize = eMean^2/(eVar-eMean)
  eProb = eMean/eVar
  numPopu = as.numeric(table(mc))
  X = rep(0,leng)
  for(i in 1:nbS){
    set.seed(seed)
    X_ = rnbinom(numPopu[i],size = eSize[i], prob = eProb[i])
    X[mc==i] = X_
  }
  val.propre  <-  round(eigen(t(trans))$values,3)
  pos         <-  which(val.propre == 1.000)
  init.tmp    <-  eigen(t(trans))$vectors[,pos]
  init        <-  init.tmp / sum(init.tmp)
  init        <-  as.numeric(init)
  return(list(x=X, mc=mc,transPr=trans,initPr=init,eParam=cbind(eMean,eVar)))
}


# leng  = 500
# trans = matrix(c(0.8, 0.2,
#          0.4, 0.6),2,2,byrow=T)
# eMean = c(100,200)
# eVar = c(1000, 1000)
# sim = simHMM(leng, trans, eMean, eVar)

