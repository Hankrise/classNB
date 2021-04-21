#Objective function
objfun = function(x){
  f = (x[2]-x[1]^2)^2 + (1-x[1])^2
  return(f)
}

#Nelder-Mead
Nelder.Mead=function(objfun, x, eps = 1.0e-6){
  t = x
  n = length(x) 
  x = matrix(0, n+1, n) 
  y = matrix(0, 1, n+1)
  x[1, ] = t
  for (i in 1:n){
    x[i+1, ] = x[1, ]
    if(x[1,i] == 0){
      x[i+1,i] = 0.00025
    }
    else{
      x[i+1,i] = 1.05*x[1,i]
    }
  }
  for (j in 1:10000){
    for (k in 1:(n+1)){
      y[k] = objfun(x[k, ])
    }
    y.order = order(y)
    y = sort(y)
    x = x[y.order, ]   
    if (norm(as.matrix(x[n+1, ] - x[1, ]),"2") < eps){
      break
    }
    m = apply(x[1:n, ], 2, mean)
    r = 2*m - x[n+1, ]
    fr = objfun(r)
    if (y[1] <= fr && fr < y[n]){
      x[n+1, ] = r
      next
    }
    else if (fr < y[1]){
      s = m + 2*(m-x[n+1, ])
      if (objfun(s) < fr){
        x[n+1, ] = s
      }
      else{
        x[n+1, ] = r
      }
      next
    }
    else if (fr < y[n+1]){
      c1 = m + (r-m)*0.5
      if (objfun(c1) < fr){
        x[n+1, ] = c1
        next
      }
    }
    else{
      c2 = m + (x[n+1, ]-m)*0.5
      if (objfun(c2) < y[n+1]){
        x[n+1, ] = c2
        next
      }
    }
    for (l in 2:(n+1)){
      x[l, ] = x[1, ] + (x[l, ]-x[1, ])*0.5
    }
  }
  output.xmin = x[1, ]
  output.fmin = objfun(x[1, ])
  output = list(xmin = output.xmin, fmin = output.fmin)
  return(output)
}

