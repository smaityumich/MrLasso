# Aggrigated estimator
Loss <- function(x, eta)
  return(min(x^2,eta^2))

mr.lasso <- function(mu, eta, accuiracy = 5e-2)
{
  by.acc = (max(mu) - min(mu)) * accuiracy
  if(by.acc > 0)
  {
    mu0 <- c(seq(min(mu),max(mu), by = accuiracy),max(mu))
    LOSS <- rep(0,length(mu0))
    for(i in 1:length(mu0))
    {
      s = 0
      for(j in 1:length(mu))
        s = s + Loss(mu0[i]-mu[j],eta)
      LOSS[i] = s
    }
    #print(which.min(LOSS))
    return(mu0[which.min(LOSS)])
  }
  else{
    return(max(mu))
  }
}


mr.lasso.i.eta <- function(obj, eta)
{
  organs = names(obj)
  i.nzero = c()
  for(organ in organs)
    i.nzero = c(i.nzero, (which(obj[[organ]]$beta[-1] != 0) + 1))
  
  i.nzero = unique(i.nzero)
  beta = obj[[organs[1]]]$beta
  v.temp = rep(0, length(organs))
  for(i.organ in 1:length(organs))
    v.temp[i.organ] = obj[[organs[i.organ]]]$beta[1]
  beta[1] = mean(v.temp)
  
  if(length(i.nzero) > 0)
  {
    locals.temp = matrix(0, length(i.nzero), length(organs))
    for(i.organ in 1:length(organs))
      locals.temp[, i.organ] = obj[[organs[i.organ]]]$beta[i.nzero]
    
    mrlasso.temp = rep(0, length(i.nzero))
    for(i in 1:length(i.nzero))
      mrlasso.temp[i] = mr.lasso(locals.temp[i, ], eta = eta)
    
    beta[i.nzero] = mrlasso.temp
  }
  
  
  ret.obj = list(beta = beta)
  for(organ in organs)
    ret.obj[[organ]] = obj[[organ]]$beta - ret.obj$beta
  
  return(ret.obj)
  
}


mr.lasso.eta <- function(beta.cv, eta)
{
  ret.obj = list()
  nfolds = length(beta.cv)
  for(i.cv in 1:nfolds)
    ret.obj[[i.cv]] = mr.lasso.i.eta(beta.cv[[i.cv]], eta = eta)
  
  return(ret.obj)
}

mr.lasso.all <- function(beta.cv, etas = c(0.1, 0.2))
{
  ret.obj = list(etas = etas, betas = list())
  for(i.eta in 1:length(etas))
  {
    eta = etas[i.eta]
    ret.obj[['betas']][[i.eta]] = mr.lasso.eta(beta.cv, eta = eta)
  }
  return(ret.obj)
}


cv.integrated.mrlasso.eta.lambda.i <- function(obj, total, test.index.majority,
                                               index.minority, lambda, index.cv = 1)
{
  organs = names(obj)
  x = as.matrix(total[, -(1:2)])
  y = total$ActArea
  
  n.organ = rep(0, length(organs))
  for(i in 2:length(organs))
  {
    organ = organs[i]
    n.organ[i] = length(which(total$Organ == organ))
  }
  
  n.organ[1] = sum(n.organ[-1])
  thresholds = lambda / sqrt(n.organ)
  
  
  for(i in 1:length(organs))
  {
    organ = organs[i]
    beta = obj[[organ]]
    i.nzero = which(beta[-1] != 0) + 1
    obj[[organ]][i.nzero] = soft.th(beta[i.nzero], thresholds[i])
    
  }
  
  errors = list()
  errors[['majority']] = list()
  for(organ in names(test.index.majority))
  {
    temp.obj = list(beta = obj$beta + obj[[organ]])
    organ.test = test.index.majority[[organ]][index.cv, ]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    error.test = y.test - predict.new(temp.obj, x.test)
    errors[['majority']][[organ]] = mean(error.test ^2 )
  }
  
  errors[['minority']] = list()
  for(organ in names(index.minority))
  {
    organ.test = index.minority[[organ]]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    error.test = y.test - predict.new(obj, x.test)
    errors[['minority']][[organ]] = mean(error.test ^2)
  }
  return(errors)
}

cv.integrated.mrlasso.eta.lambda <- function(obj, total, test.index.majority,
                                             index.minority, lambda)
{
  organs = names(test.index.majority)
  nfold = nrow(test.index.majority[[organs[1]]])
  
  
  cv.mse = list()
  
  for(i.cv in 1:nfold)
  {
    cv.mse[[i.cv]] =  cv.integrated.mrlasso.eta.lambda.i(obj[[i.cv]], total, test.index.majority,
                                                         index.minority, lambda, i.cv)
  }
  
  return.mse = list()
  
  for(cat in c('majority', 'minority'))
  {
    organs = names(cv.mse[[1]][[cat]])
    mse = matrix(0, length(organs), nfold)
    for(i.cv in 1:nfold)
    {
      mse[, i.cv] = unlist(cv.mse[[i.cv]][[cat]])
    }
    
    row.names(mse) = organs
    return.mse[[cat]] = mse
  }
  return(return.mse)
}
