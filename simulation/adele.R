get.adele.dense <- function(obj)
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
  
  beta[i.nzeros] = apply(beta.local, 1, mean)
  return(list(beta = beta))
}


get.adele.cv <- function(local.obj, adele.dense, data.valid, lambda = 1)
{
  t.local = lambda
  t.global = lambda / sqrt(data.valid$pars$m)
  adele.sparse = soft.th(adele.dense$beta, t.global)
  
  delta.local = list()
  for(k in 1:length(local.obj))
    delta.local[[k]] = soft.th(local.obj[[k]]$beta - adele.dense$beta, t.local)
  
  error = 0
  for(k in 1:length(local.obj))
  {
    sp.beta = adele.sparse + delta.local[[k]]
    error.k = data.valid$y[[k]] -  predict.new(list(beta = sp.beta), data.valid$x[[k]])
    error = error + mean(error.k ^2)
  }
  return(error/length(local.obj))
}




