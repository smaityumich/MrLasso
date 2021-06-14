setwd('~/MRLASSO/simulation')
source('data.R')
source('lasso.R')
source('utils.R')
source('adele.R')
source('mrlasso.R')
library(rjson)



instance.i <- function(n = 100, m = 5, d = 2000, s = 5, noise = 0.05)
{
  data.train = get.data(m, d, n, s, noise)
  # data.valid = get.data(m, d, 100, s, noise)
  
  global.obj = global.lasso(data.train)
  l2.global = l2.norm(global.obj$beta - data.train$beta.adele)
  
  local.obj = local.lasso(data.train)
  
  ## ADELE
  
  adele.dense = get.adele.dense(local.obj)
  l2.adele.dense = l2.norm(adele.dense$beta - data.train$beta.adele)
  
  t.global = l.infty.norm(adele.dense$beta - data.train$beta.adele)
  
 
  adele.sparse = soft.th(adele.dense$beta, t.global)
  
  l2.adele.sparse = l2.norm(adele.sparse - data.train$beta.adele)
  
  ret.obj = list(m = m, n = n, d = d, noise = noise, 
                 l2.global = l2.global, l2.adele.dense = l2.adele.dense,
                 l2.adele.sparse = l2.adele.sparse)
  
  ## MRLASSO
  
  etas = c(2, 4, 8, 16, 32, 64)
  
  for(i in 1:length(etas))
  {
    eta = etas[i]
    beta = data.train$beta
    
    beta.mrlasso = data.train$beta.mrlasso
    ml.eta <- function(x)
    {
      return(mr.lasso(x, eta = eta))
    }
    beta.mrlasso[2:(d+1)] = apply(beta, 1, ml.eta)
    
    mrlasso.dense = get.mrlasso.dense(local.obj, eta = eta)
    l2.mrlasso.dense = l2.norm(mrlasso.dense$beta - beta.mrlasso)
    
    t.global = l.infty.norm(mrlasso.dense$beta - beta.mrlasso)
    mrlasso.sparse = soft.th(mrlasso.dense$beta, t.global)
    l2.mrlasso.sparse = l2.norm(mrlasso.sparse - beta.mrlasso)
    
    ret.obj[[paste('l2.mrlasso.dense.', eta, sep = '')]] = l2.mrlasso.dense
    ret.obj[[paste('l2.mrlasso.sparse.', eta, sep = '')]] = l2.mrlasso.sparse
  }
  
  
  
  

 
  
  
  
  return(ret.obj)
}



instance <- function(n = 100, m = 5, d = 2000, s = 1, noise = 0.25, iter = 50)
{
  l2 = matrix(0, iter, 5)
  l2 = data.frame(l2)
  colnames(l2) <- c('l2.global', 'l2.adele.dense',
                    'l2.adele.sparse', 'l2.mrlasso.dense',
                    'l2.mrlasso.sparse', 'm', 'n', 'sigma')
  for(i in 1:iter)
  {
    e = instance.i(n = n, m = m, d = d, s = s, noise = noise)
    l2[i, 1] = e$l2.global
    l2[i, 2] = e$l2.adele.dense
    l2[i, 3] = e$l2.adele.sparse
    l2[i, 4] = e$l2.mrlasso.dense
    l2[i, 5] = e$l2.mrlasso.sparse
    l2[i, 6] = e$m
    l2[i, 7] = e$n
    l2[i, 8] = e$noise
  }
  return(l2)
}


v.sigma = c(rep(0.1, 10), rep(0.25, 10))
v.n = c(100 * 2^(0:4), rep(100, 5))
v.n = c(v.n, v.n)
v.m = c(rep(5, 5) ,2^(1:5))
v.m = c(v.m, v.m)

pars = data.frame(sigma = v.sigma, m = v.m, n = v.n)





args = commandArgs(trailingOnly=TRUE)
i = as.integer(args[1])


index = i %/% 50 + 1
sig = pars[index, 1]
m = pars[index, 2]
n = pars[index, 3]

set.seed(i)

l2 = instance.i(n = n, m = m, noise = sig, s = 5)

write(rjson::toJSON(l2), file = paste('l2/l2_', i, '.json', sep = ''))



