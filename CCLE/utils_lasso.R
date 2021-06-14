library(glmnet)


lasso <- function(x, y, lambda = NULL, nlambda = 10, nfolds = 7)
{
  #lambda = NULL; nlambda = 10; nfolds = 5
  if(is.null(lambda))
  {
    cv = cv.glmnet(x, y, nfolds = nfolds,
                   trace.it = T, standardize = T, nlambda = nlambda)
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
    tau.i = mean(error.i^2) + obj.i$lambda * sum(abs(obj.i$beta))
    obj$beta[i+1] = obj$beta[i+1] + mean(error.i * error.y)/tau.i
  }
  return(obj)
}




global.lasso.cv.i <- function(total, majority.total.index, test.index.majority, index.minority, index.cv = 1)
{
  test.index.total = c()
  for(organ in test.index.majority)
    test.index.total = c(test.index.total, organ[index.cv, ])
  train.index.total = setdiff(majority.total.index, test.index.total)
  
  y = total$ActArea
  x = as.matrix(total[, -(1:2)])
  y.train = y[train.index.total]
  x.train = x[train.index.total, ]
  global = lasso(x.train, y.train)
  
  
  
  errors = list()
  
  errors[['majority']] = list()

  for(organ in names(test.index.majority))
  {
    organ.test =test.index.majority[[organ]][index.cv, ]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    error.test = y.test - predict.new(global, x.test)
    errors[['majority']][[organ]] = error.test
  }
  
  errors[['minority']] = list()
  for(organ in names(index.minority))
  {
    organ.test = index.minority[[organ]]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    error.test = y.test - predict.new(global, x.test)
    errors[['minority']][[organ]] = error.test
  }
  return(errors)
}


global.mean.cv.i <- function(total, majority.total.index, test.index.majority, index.minority, index.cv = 1)
{
  
  
  y = total$ActArea
  x = as.matrix(total[, -(1:2)])
  
  local.means = list()
  for(organ in names(test.index.majority))
  {
    index.organ = which(total$Organ == organ)
    index.test = test.index.majority[[organ]][index.cv, ]
    index.train = setdiff(index.organ, index.test)
    y.train = y[index.train]
  
    local.means[[organ]] = mean(y.train)
  }
  
  
  global.mean = mean(as.numeric(unlist(local.means)))
  
  errors = list()
  errors[['majority']] = list()
  for(organ in names(test.index.majority))
  {
    organ.test =test.index.majority[[organ]][index.cv, ]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    errors[['majority']][[organ]] = y.test - global.mean
  }
  
  errors[['minority']] = list()
  for(organ in names(index.minority))
  {
    organ.test = index.minority[[organ]]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    errors[['minority']][[organ]] = y.test - global.mean
  }
  return(errors)
}






global.lasso <- function(total, majority.total.index, test.index.majority, index.minority)
{
  n = names(test.index.majority)
  if(length(n) > 0)
  {
    nfold = nrow(test.index.majority[[n[1]]])
    cv = list()
    for(i.cv in 1:nfold)
      cv[[i.cv]] = global.lasso.cv.i(total, majority.total.index, test.index.majority, index.minority, i.cv)
    
    
    return.mse = list()
    
    for(cat in c('majority', 'minority'))
    {
      organs = names(cv[[1]][[cat]])
      mse = matrix(0, length(organs), nfold)
      for(i.organ in 1:length(organs))
      {
        organ = organs[i.organ]
        #errs = c()
        for(i.cv in 1:nfold)
        {
          mse[i.organ, i.cv] = mean(cv[[i.cv]][[cat]][[organ]] ^2 )
        }
        #mse[i.organ] = mean(errs ^ 2)
      }
      row.names(mse) = organs
      return.mse[[cat]] = mse
    }
    return(return.mse)
  }
  else
  {
    return(NULL)
  }
}


global.mean <- function(total, majority.total.index, test.index.majority, index.minority)
{
  n = names(test.index.majority)
  if(length(n) > 0)
  {
    nfold = nrow(test.index.majority[[n[1]]])
    cv = list()
    for(i.cv in 1:nfold)
      cv[[i.cv]] = global.mean.cv.i(total, majority.total.index, test.index.majority, index.minority, i.cv)
    
    
    return.mse = list()
    
    for(cat in c('majority', 'minority'))
    {
      organs = names(cv[[1]][[cat]])
      mse = matrix(0, length(organs), nfold)
      for(i.organ in 1:length(organs))
      {
        organ = organs[i.organ]
        #errs = c()
        for(i.cv in 1:nfold)
        {
          mse[i.organ, i.cv] = mean(cv[[i.cv]][[cat]][[organ]] ^2 )
        }
        #mse[i.organ] = mean(errs ^ 2)
      }
      row.names(mse) = organs
      return.mse[[cat]] = mse
    }
    return(return.mse)
  }
  else
  {
    return(NULL)
  }
}



local.lasso.cv.i <- function(total, test.index.majority, index.cv = 1)
{
  local.pars = list()
  x = as.matrix(total[, -(1:2)])
  y = total$ActArea
  
  for(organ in names(test.index.majority))
  {
    index.organ = which(total$Organ == organ)
    index.test = test.index.majority[[organ]][index.cv, ]
    index.train = setdiff(index.organ, index.test)
    y.train = y[index.train]
    x.train = x[index.train, ]
    obj = lasso(x.train, y.train)
    obj = lasso.debias(x.train, y.train, obj, nfolds = 7)
    local.pars[[organ]] = obj
  }
  return(local.pars)
}


local.lasso <- function(total, test.index.majority)
{
  n = names(test.index.majority)
  if(length(n) > 0)
  {
    nfold = nrow(test.index.majority[[n[1]]])
    local.parameter.cv = list()
    for(i.cv in 1:nfold)
      local.parameter.cv[[i.cv]] = local.lasso.cv.i(total, test.index.majority, index.cv = i.cv)
    
    return(local.parameter.cv)
  }
  else
  {
    return(NULL)
  }
}

