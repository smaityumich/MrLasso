library(glmnet)


adele.i <- function(obj, total, lambda)
{
  
  organs = names(obj)
  beta = obj[[organs[1]]]$beta
  for(i in 2:length(organs))
  {
    beta = beta + obj[[organs[i]]]$beta
  }
  beta = beta / length(organs)
  ret.obj = list(beta = beta)
  for(organ in organs)
    ret.obj[[organ]] = obj[[organ]]$beta - ret.obj$beta
  
  
  n.organ = rep(0, length(organs))
  for(i in 1:length(organs))
  {
    organ = organs[i]
    n.organ[i] = length(which(total$Organ == organ))
  }
  
  thresholds = c(sum(n.organ), n.organ)
  thresholds = lambda / sqrt(thresholds)
  
  beta = ret.obj$beta
  i.nzero = which(beta[-1] != 0) + 1
  
  ret.obj$beta[i.nzero] = soft.th(beta[i.nzero], thresholds[1])
  
  for(i in 1:length(organs))
  {
    organ = organs[i]
    beta = ret.obj[[organ]]
    i.nzero = which(beta[-1] != 0) + 1
    ret.obj[[organ]][i.nzero] = soft.th(beta[i.nzero], thresholds[i+1])
    
  }
  
  
  
  
  return(ret.obj)
}


mse.adele.i <- function(obj, total, test.index.majority, index.minority, lambda, index.cv = 1)
{
  x = as.matrix(total[, -(1:2)])
  y = total$ActArea
  ret.obj = adele.i(obj, total, lambda)
  
  errors = list()
  errors[['majority']] = list()
  for(organ in names(test.index.majority))
  {
    temp.obj = list(beta = ret.obj$beta + ret.obj[[organ]])
    organ.test = test.index.majority[[organ]][index.cv, ]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    error.test = y.test - predict.new(temp.obj, x.test)
    errors[['majority']][[organ]] = error.test
  }
  
  errors[['minority']] = list()
  for(organ in names(index.minority))
  {
    organ.test = index.minority[[organ]]
    y.test = y[organ.test]
    x.test = x[organ.test, ]
    error.test = y.test - predict.new(ret.obj, x.test)
    errors[['minority']][[organ]] = error.test
  }
  return(errors)
}

mse.adele <- function(beta.cv, total, test.index.majority, index.minority, lambda)
{
  organs = names(test.index.majority)
  
  nfold = nrow(test.index.majority[[organs[1]]])
  
  cv = list()
  for(i.cv in 1:nfold)
  {
    cv[[i.cv]] = mse.adele.i(beta.cv[[i.cv]], total, test.index.majority, index.minority, lambda, index.cv = i.cv)
  }
  
  
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



