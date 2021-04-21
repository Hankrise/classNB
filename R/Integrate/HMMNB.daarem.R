HMMNB.daarem <- function(par, fixptfn, objfn, ..., control=list()) {

  control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0, kappa=25, alpha=1.2)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)

  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa

  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))

  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)

  xold <- par
  xnew <- fixptfn(xold, ...)
  obj_funvals[1] <- objfn(xold, ...)
  obj_funvals[2] <- objfn(xnew, ...)
  likchg <- obj_funvals[2] - obj_funvals[1]
  obj.evals <- 2

  fold <- xnew - xold
  k <- 1
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- TRUE
  #num.em <- 0  ## number of EM fallbacks
  ell.star <- obj_funvals[2]
  while(k < maxiter) {
    count <- count + 1

    fnew <- fixptfn(xnew, ...) - xnew
    ss.resids <- sqrt(crossprod(fnew))

    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold
    if(ss.resids < tol & count==nlag) break


    np <- count
    Ftmp <- matrix(Fdiff[,1:np], nrow=length(fnew), ncol=np)
    Xtmp <- matrix(Xdiff[,1:np], nrow=length(fnew), ncol=np)  ## is this matrix function needed?

    tmp <- svd(Ftmp)
    dvec <- tmp$d
    uy <- crossprod(tmp$u, fnew)
    uy.sq <- uy*uy

    ### Still need to compute Ftf
    Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew))^2))
    tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr

    dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
    gamma_vec <- tmp$v%*%dd

    if(class(gamma_vec) != "try-error"){

      xbar <- xnew - drop(Xtmp%*%gamma_vec)
      fbar <- fnew - drop(Ftmp%*%gamma_vec)

      x.propose <- xbar + fbar
      new.objective.val <- try(objfn(x.propose, ...), silent=TRUE)
      obj.evals <- obj.evals + 1

      if(class(new.objective.val) != "try-error" & !is.na(obj_funvals[k+1]) &
         !is.nan(new.objective.val)) {
        if(new.objective.val >= obj_funvals[k+1] - mon.tol) {
          ## Increase delta
          obj_funvals[k+2] <- new.objective.val
          fold <- fnew
          xold <- xnew

          xnew <- x.propose
          shrink.count <- shrink.count + 1
        } else {
          ## Keep delta the same
          fold <- fnew
          xold <- xnew

          xnew <- fold + xold
          obj_funvals[k+2] <- objfn(xnew, ...)
          obj.evals <- obj.evals + 1
          #num.em <- num.em + 1
        }
      } else {
        ## Keep delta the same
        fold <- fnew
        xold <- xnew

        xnew <- fold + xold
        obj_funvals[k+2] <- objfn(xnew, ...)
        obj.evals <- obj.evals + 1
        count <- 0
        #num.em <- num.em + 1
      }
    } else {
      ## Keep delta the same
      fold <- fnew
      xold <- xnew

      xnew <- fold + xold
      obj_funvals[k+2] <- objfn(xnew, ...)
      obj.evals <- obj.evals + 1
      count <- 0
      #num.em <- num.em + 1
    }
    if(count==nlag) {
      count <- 0
      ## restart count
      ## make comparison here l.star vs. obj_funvals[k+2]
      if(obj_funvals[k+2] < ell.star - cycl.mon.tol) {
        ## Decrease delta
        shrink.count <- max(shrink.count - nlag, -2*kappa)
      }
      ell.star <- obj_funvals[k+2]
    }

    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
    k <- k+1
  }
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]
  value.obj <- objfn(xnew, ...)
  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, convergence=conv, objfn.track=obj_funvals))
}
HMMNBLogLik <- function(para.hat, X, nbS) {
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
  # nbT = ncol(Fpr); nbS = nrow(Fpr) # backward
  # postPr = matrix(NA, nbS, nbT)
  # Gpr = matrix(NA, nbS, nbT)
  # postPr = matrix(NA, nbS, nbT)
  # for(tt in (nbT-1):1){
  #   Gpr[,tt+1]  = Fpr[,tt] %*% transPr
  #   postPr[,tt] = Fpr[,tt] *(transPr %*% (postPr[,tt+1]/Gpr[,tt+1]))#sapply(1:nbS, function(l) sum(trans[l,] * postPr[,tt+1] / Gpr[,tt+1]))#
  # }
  # loglik.a=sum(postPr[,1]*log(initPr))
  # loglik.b=
  # loglik.c=sum(postPr*Pmat)
  # loglik  = loglik.a+loglik.b+loglik.c
  return(logLik)
}
DampingFind <- function(uy.sq, dvec, aa, kappa, sk, Ftf, lambda.start=NULL,
                        r.start=NULL, maxit=10) {
  ## uy.sq is a vector whose jth component (U^t*y)^2
  ## d.sq is a vector whose jth component is d_j^2.
  ## (aa, kappa, sk) - parameters in determining delta_k
  ## Ftf <- F_{k}^T f_{k}
  ##
  ## Output: value of penalty term such that (approximately)
  ## .       || beta.ridge ||^2/||beta.ols||^2 = delta_{k}^{2}, with 0 < target < 1.

  if(is.null(lambda.start) | is.na(lambda.start)) {
    lambda.start <- 100
  }
  if(is.null(r.start) | is.na(r.start)) {
    r.start <- 0
  }
  if(sum(dvec==0.0) > 0) {
    ind <- dvec > 0
    dvec <- dvec[ind]
    uy.sq <- uy.sq[ind]
    ## what about if dvec has length 1?
  }
  pow <- kappa - sk
  target <- exp(-0.5*log1p(aa^pow))  ### This is sqrt(delta_k)
  #d.sq <- pmax(dvec*dvec, 1e-7)  ## in case the least-squares problem is very poorly conditioned
  d.sq <- dvec*dvec
  betahat.ls <- uy.sq/d.sq
  betahat.ls.norm <- sqrt(sum(betahat.ls))
  vk <- target*betahat.ls.norm

  if(vk == 0) {
    ## the norm of the betas is zero, so the value of lambda shouldn't matter
    return(list(lambda=lambda.start, rr=r.start))
  }
  ### Initialize lambda and lower and upper bounds
  lambda <- lambda.start - r.start/vk
  LL <- (betahat.ls.norm*(betahat.ls.norm - vk))/sum(uy.sq/(d.sq*d.sq))
  UU <- Ftf/vk

  ## Compute lstop and ustop
  pow.low <- pow + 0.5
  pow.up <- pow - 0.5
  l.stop <- exp(-0.5*log1p(aa^pow.low))
  u.stop <- exp(-0.5*log1p(aa^pow.up))

  ### Start iterations
  for(k in 1:maxit) {
    ### Check if lambda is inside lower and upper bounds
    if(lambda <= LL | lambda >= UU) {
      lambda <- max(.0001*UU, sqrt(LL*UU))
    }
    ### Evaluate ||s(lambda)|| and \phi(lambda)/phi'(lambda)

    d.lambda <- (dvec/(d.sq + lambda))^2
    d.prime <- d.lambda/(d.sq + lambda)
    s.norm <- sqrt(sum(uy.sq*d.lambda))
    phi.val <- s.norm - vk
    phi.der <- (-1)*sum(uy.sq*d.prime)/s.norm
    phi.ratio <- phi.val/phi.der

    d.u <- (dvec/(d.sq + LL))^2
    s.up <- sqrt(sum(uy.sq*d.u))
    ### Initial Lower bound is not correct

    ### Check convergence
    if(s.norm <= u.stop*betahat.ls.norm & s.norm >= l.stop*betahat.ls.norm) {
      break
    }

    ### If not converged, update lower and upper bounds
    UU <- ifelse(phi.val >= 0, UU, lambda)
    LL <- max(LL, lambda - phi.ratio)

    ### Now update lambda
    lambda <- lambda - (s.norm*phi.ratio)/vk
    bb <- (s.norm*phi.ratio)/vk
  }
  ans <- list(lambda=lambda, rr=s.norm*phi.ratio)
  return(ans)
}
