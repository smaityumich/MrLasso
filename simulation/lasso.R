library(glmnet)

lasso <- function(x, y, lambda = NULL, nlambda = 10, nfolds = 5)
{
  
  if(is.null(lambda))
  {
    cv = cv.glmnet(x, y, nfolds = nfolds,
                   trace.it = F, standardize = T, nlambda = nlambda)
    ind.best = which(cv$lambda == cv$lambda.min)
    
    coeff = cv$glmnet.fit$beta[, ind.best]
    lambda.best = cv$lambda.min
    coeff.best = coef.glmnet(cv, s = 'lambda.min')
    intercept.best = cv$glmnet.fit$a0[ind.best]
  }
  
  else
  {
    lambda.best = lambda
    cv = glmnet(x, y, family = 'gaussian', lambda = lambda.best)
    coeff.best = coef.glmnet(cv)
  }
  return(list(lambda = lambda.best, beta = coeff.best))
}


predict.new <- function(obj, newx)
{
  return(newx %*% obj$beta[-1] + obj$beta[1])
}




lasso.debias <- function(x, y, obj, nlambda = 10, nfolds = 3)
{
  beta = obj$beta[-1]
  error.y = y - predict.new(obj, x)
  
  nzeros = which(beta != 0)
  for(i in nzeros)
  {
    y.i = x[, i]
    x.i = x[, -i]
    obj.i = lasso(x.i, y.i, nlambda = nlambda, nfolds = nfolds)
    error.i = y.i - predict.new(obj.i, x.i)
    tau.i.sq = mean(error.i^2) + obj.i$lambda * sum(abs(obj.i$beta))
    obj$beta[i+1] = obj$beta[i+1] + mean(error.i * error.y)/tau.i.sq
  }
  
  zeros = which(beta == 0)
  zeros.sample = sample(zeros, 20)
  replace.debias = rep(0, 20)


  for(j in 1:20)
  {
    i = zeros.sample[j]
    y.i = x[, i]
    x.i = x[, -i]
    obj.i = lasso(x.i, y.i, nlambda = nlambda, nfolds = nfolds)
    error.i = y.i - predict.new(obj.i, x.i)
    tau.i.sq = mean(error.i^2) + obj.i$lambda * sum(abs(obj.i$beta))
    obj$beta[i+1] = obj$beta[i+1] + mean(error.i * error.y)/tau.i.sq
    replace.debias[j] =  obj$beta[i+1]
  }

  zeros.not.in.sample = setdiff(zeros, zeros.sample) + 1
  obj$beta[zeros.not.in.sample] = sample(replace.debias, length(zeros.not.in.sample), replace = T)
  
  
  return(obj)
}


global.lasso <- function(data)
{
  print('Running global lasso ...')
  x.combined = data$x[[1]]
  for(k in 2:data$pars$m)
    x.combined = rbind(x.combined, data$x[[k]])
  
  y.combined = data$y[[1]]
  for(k in 2:data$pars$m)
    y.combined = rbind(y.combined, data$y[[k]])
  
  obj = lasso(x.combined, y.combined, nfolds = 10)
  return(obj)
}


local.lasso <- function(data)
{
  obj = list()
  for(k in 1:data$pars$m)
  {
    print(paste('Running local lasso for data', k, '...'))
    temp = lasso(data$x[[k]], data$y[[k]], nlambda = 10, nfolds = 5)
    obj[[k]] = lasso.debias(data$x[[k]], data$y[[k]], obj = temp,  nlambda = 10, nfolds = 5)
  }
  return(obj)  
}



