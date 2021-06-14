

local.lasso.final <- function(total, organs, lambda.adele.best, lambda.best.mrlasso, eta.mrlasso.best)
{
  local.pars = list()
  x = as.matrix(total[, -(1:2)])
  y = total$ActArea
  
  local.pars[['beta']] = list()
  debiased = list()
  
  for(organ in organs)
  {
    index.organ = which(total$Organ == organ)

    y.organ = y[index.organ]
    x.organ = x[index.organ, ]
    obj = lasso(x.organ, y.organ)
    local.pars[['beta']][[organ]] = obj
    obj = lasso.debias(x.organ, y.organ, obj, nfolds = 7)
    debiased[[organ]] = obj
  }
  
  n.majority = sum(total$Organ %in% organs)
  
  #ADELE
  temp.obj = adele.i(debiased, total, lambda = lambda.adele.best/sqrt(n.majority))
  local.pars[['adele']] = temp.obj$beta
  
  #MRLASSO
  temp.obj = mr.lasso.i.eta(debiased, eta = eta.mrlasso.best)[['beta']]
  i.nzero = which(temp.obj[-1] != 0) + 1
  temp.obj[i.nzero] = soft.th(temp.obj[i.nzero], lambda.best.mrlasso/sqrt(n.majority))
  
  local.pars[['mrlasso']] = temp.obj
  
  #Global.LASSO
  index.total = which(total$Organ %in% organs)
  y.organ = y[index.total]
  x.organ = x[index.total, ]
  obj = lasso(x.organ, y.organ)
  local.pars[['global.lasso']] = obj
  
  
  return(local.pars)
}


