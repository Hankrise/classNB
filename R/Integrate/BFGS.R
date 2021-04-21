#Objective function
objfun = function(x){
  f = (x[2]-x[1]^2)^2 + (1-x[1])^2
  return(f)
}

#Gradient of the objective function
gradfun=function(x){
  g = as.matrix(c(4.0*(x[1]^3-x[1]*x[2]) + 2*x[1]-2, 2.0*(x[2]-x[1]^2)))
  return(g)
}

BFGS = function(objfun, gradfun, x_initial, eps = 1.0e-10){
  ite = 0


  x = x_initial
  n = length(x)
  H = diag(n)
  while(ite < 10000){
    d = -H%*%gradfun(x)
    lambda = argst(objfun, x, d)
    x_last = x
    x = x + lambda*d
    delta_g = gradfun(x) - gradfun(x_last)
    delta_x = x - x_last
    # print(delta_g)
    # print(delta_x)
    if (norm(delta_x)<eps || norm(delta_g)<eps){
      break
    }
    H = (diag(n)-(delta_x%*%t(delta_g))/as.numeric(t(delta_x)%*%delta_g)) %*% H %*% (diag(n)-(delta_g%*%t(delta_x))/as.numeric(t(delta_x)%*%delta_g)) + (delta_x%*%t(delta_x))/as.numeric(t(delta_x)%*%delta_g)
    ite = ite + 1
  }
  output.xmin = x[ ,1]
  output.fmin = objfun(x)
  output = list(xmin = output.xmin, fmin = output.fmin)
  return(output)
}

#Line search with Golden-section
argst=function(objfun, x, d, t0 = 0, step = 0.0001, eps = 1e-06){
  t1 = t0 + step
  ft0 = objfun(x + t0*d)
  ft1 = objfun(x + t1*d)
  if(ft1 <= ft0){
    step = 2*step
    t2 = t1 + step
    ft2 = objfun(x + t2*d)
    while(ft1 > ft2){
      t1 = t2
      step = 2*step
      t2 = t1 + step
      ft1 = objfun(x + t1*d)
      ft2 = objfun(x + t2*d)
    }
  }else{
    step = step/2
    t2 = t1
    t1 = t2 - step
    ft1 = objfun(x + t1*d)
    while(ft1 > ft0){
      step = step/2
      t2 = t1
      t1 = t2 - step
      ft1 = objfun(x + t1*d)
    }
  }
  a = 0
  b = t2
  a1 = a + 0.382*(b - a)
  a2 = a + 0.618*(b - a)
  f1 = objfun(x + a1*d)
  f2 = objfun(x + a2*d)
  while(abs(b - a) >= eps){
    if(f1 < f2){
      b = a2
      a2 = a1
      f2 = f1
      a1 = a + 0.382*(b - a)
      f1 = objfun(x + a1*d)
    }else{
      a = a1
      a1 = a2
      f1 = f2
      a2 = a + 0.618*(b - a)
      f2 = objfun(x + a2*d)
    }
  }
  alpha = 0.5*(a + b)
  return(alpha)
}

