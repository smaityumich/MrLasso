mr.lasso <- function(x, eta, tol = 1e-9, max.iter = 1000)
{
  # if(sum(x == 0) >= length(x)/2)
  #   return(0)
  if(T)
  {
    mu0 <- mean(x) #+ rnorm(1) * 0.05
    error = 1
    iter = 1
    while(error>tol & iter < max.iter)
    {
      diff = abs(x-mu0)
      y = x[which(diff<=eta + abs(mu0))]
      mu1 = mean(y) #+ rnorm(1) * 0.05 / (iter ^ 2)
      error = abs(mu0-mu1) #
      mu0=mu1
      if(is.na(mu0))
        mu0 = mean(x) +  rnorm(1) * 0.05
      iter = iter+1
    }
    return(mu0)
  }
  
}



get.mrlasso.dense <- function(obj, eta = 0.4)
{
  d = obj[[1]]$beta@Dim[1]
  beta = Matrix::Matrix(0, d, 1)
  
  i.nzeros = c()
  for(o in obj)
  {
    i.nzeros = c(i.nzeros, which(o$beta != 0))
  }
  i.nzeros = sort(unique(i.nzeros))
  
  beta.local = matrix(0, length(i.nzeros), length(obj))
  for(k in 1:length(obj))
    beta.local[, k] = obj[[k]]$beta[i.nzeros]
  
  ml.eta <- function(x)
  {
    return(mr.lasso(x, eta = eta))
  }
  
  beta[i.nzeros] = apply(beta.local, 1, ml.eta)
  return(list(beta = beta))
}


get.mrlasso.cv <- function(local.obj, mrlasso.dense, data.valid, lambda = 1)
{
  t.local = lambda
  t.global = lambda / sqrt(data.valid$pars$m)
  mrlasso.sparse = soft.th(mrlasso.dense$beta, t.global)
  
  delta.local = list()
  for(k in 1:length(local.obj))
    delta.local[[k]] = soft.th(local.obj[[k]]$beta - mrlasso.dense$beta, t.local)
  
  error = 0
  for(k in 1:length(local.obj))
  {
    sp.beta = mrlasso.sparse + delta.local[[k]]
    error.k = data.valid$y[[k]] -  predict.new(list(beta = sp.beta), data.valid$x[[k]])
    error = error + mean(error.k ^2)
  }
  return(error/length(local.obj))
}





eval.mrlasso <- function(eta, local.obj, beta.mrlasso)
{
  mrlasso.dense = get.mrlasso.dense(local.obj, eta = eta)
  
  l2.mrlasso.dense = l2.norm(mrlasso.dense$beta - beta.mrlasso)
  
  t.global = l.infty.norm(mrlasso.dense$beta - beta.mrlasso)
  mrlasso.sparse = soft.th(mrlasso.dense$beta, t.global)
  l2.mrlasso.sparse = l2.norm(mrlasso.sparse - beta.mrlasso)
  return(c(l2.mrlasso.dense, l2.mrlasso.sparse))
}

get.best.eta <- function(local.obj, beta.mrlasso)
{
  start = 0
  
  etas = seq(1, 10, by = 0.5)
  l2s = rep(0, length(etas))
  for(i in 1:length(etas))
     l2s[i] = eval.mrlasso(etas[i], local.obj = local.obj,
                          beta.mrlasso = beta.mrlasso)[2]
  eta = as.numeric(quantile(etas[l2s == min(l2s)], probs = 0.2))
  
  etas = seq(-5, 5, by = 1)/20 + eta
  l2s = rep(0, 10)
  for(i in 1:10)
    l2s[i] = eval.mrlasso(etas[i], local.obj = local.obj,
                          beta.mrlasso = beta.mrlasso)[2]
  eta = min(etas[l2s == min(l2s)])
  return(eta)
  
}

