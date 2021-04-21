mixNB.Sim <- function(nbT,weight,mean,variance,seed=2019)
{
  if(FALSE){
    nbT=1000
    weight=c(0.2,0.6,0.2)
    mean=c(100,200,300)
    variance=rep(1000,3)
    seed=2019

    sim=mixNB.Sim(nbT,weight,mean,variance,seed)
    sim
    x=sim$x
    print(sim,"Y")
    hist(sim)
    plot(sim)
  }
  if(!is.numeric(c(nbT,weight,mean,variance,seed))) {
    stop("Parameters have to be numeric")
  } else if(!all(c(nbT,weight,mean,variance,seed)>0)){
    stop("Parameters have to be positive")
  } else if(length(weight)!=length(mean)|length(weight)!=length(variance)){
    stop("Unmached length of Weight,Mean or Varience")
  } else {

    set.seed(seed)
    nbS=length(weight)
    eSize_set = mean^2/(variance-mean)
    eProb_set = mean/variance
    Z = sample(1:nbS, nbT, prob=weight, replace=TRUE)
    numPopu = table(Z)
    X_ = list()
    X = rep(0,nbT)
    for(i in 1:nbS){
      set.seed(seed)
      X_[[i]] = rnbinom(numPopu[i],size = eSize_set[i], prob = eProb_set[i])
      X[Z==i] = X_[[i]]
    }
    # hist(X)
    obj <-list(nbT=nbT,
               weight=weight,
               mean=mean,
               variance=variance,
               seed=seed,
               x=X,
               idx=Z)
    class(obj)<-"mixNB.Sim"
    return(obj)
  }
}
