## set project path in local machine 
project.path = "~/projects/MrLasso/"

setwd(paste(project.path, "simulation/", sep =""))
source('data.R')
source('lasso.R')
source('utils.R')
source('adele.R')
source('mrlasso.R')
library(rjson)


n = 100; m = 8; d = 2000; s = 4; s.extra = 15;
noise = 0.25; snr = 4;
outlier.dist = 10; seed = 0;



instance.i <- function(n = 100, m = 5, d = 2000, s = 5, s.extra = 15,
                       noise = 0.05, snr = 4,
                       outlier.dist = 10, seed = 0)
{
  data.train = get.data(m = m, d = d, n = n, snr = snr, 
                        noise = noise, s = s, s.extra = s.extra,
                        outlier.dist = outlier.dist, seed = seed)
  
  data.test = get.data.test(m = 2000, d = d, n = 10, snr = snr, 
                        noise = 0, s = s, s.extra = 0,
                        outlier.dist = outlier.dist, seed = seed+10000)
  # data.valid = get.data(m, d, 100, s, noise)

  
  
  global.obj = global.lasso(data.train)
  l2.global = l2.norm(global.obj$beta - data.train$beta.adele)
  # error.global = eval.central.beta(data.test, global.obj$beta)
  
  error.global = eval.central.beta.th(global.obj$beta, snr = snr, s = s, noise = noise)

  
  
  local.obj = local.lasso(data.train)
  
  ## ADELE
  
  adele.dense = get.adele.dense(local.obj)
  l2.adele.dense = l2.norm(adele.dense$beta - data.train$beta.adele)
  
  t.global = l.infty.norm(adele.dense$beta - data.train$beta.adele)
  
 
  adele.sparse = soft.th(adele.dense$beta, t.global)
 
  l2.adele.sparse = l2.norm(adele.sparse - data.train$beta.adele)
  
  # error.adele.dense = eval.central.beta(data.test, adele.dense$beta)
  # error.adele.sparse = eval.central.beta(data.test, adele.sparse)
  # 
  error.adele.dense = eval.central.beta.th(adele.dense$beta, snr = snr, s = s, noise = noise)
  error.adele.sparse = eval.central.beta.th(adele.sparse, snr = snr, s = s, noise = noise)
  
  print("adele....")
  print("true:")
  print(as.numeric(data.train$beta.adele)[1:25])
  print("estimate (sparse):")
  print(as.numeric(adele.sparse)[1:25])
  print("l2 error:")
  print(l2.adele.sparse)
  print("test data error:")
  print(error.adele.sparse)
  
  
  
  ## MRLASSO
  beta.mrlasso = data.train$beta.mrlasso
  eta = get.best.eta(local.obj, beta.mrlasso)
  
  mrlasso.dense = get.mrlasso.dense(local.obj, eta = eta)
  
  l2.mrlasso.dense = l2.norm(mrlasso.dense$beta - beta.mrlasso)
  
  t.global = l.infty.norm(mrlasso.dense$beta - beta.mrlasso)
  mrlasso.sparse = soft.th(mrlasso.dense$beta, t.global)
  l2.mrlasso.sparse = l2.norm(mrlasso.sparse - beta.mrlasso)
  
  # error.mrlasso.dense = eval.central.beta(data.test, mrlasso.dense$beta)
  # error.mrlasso.sparse = eval.central.beta(data.test, mrlasso.sparse)
    
  error.mrlasso.dense = eval.central.beta.th( mrlasso.dense$beta, snr = snr, s = s, noise = noise)
  error.mrlasso.sparse = eval.central.beta.th(mrlasso.sparse, snr = snr, s = s, noise = noise)
  
  print("MrLasso....")
  print("true:")
  print(as.numeric(data.train$beta.mrlasso)[1:25])
  print("estimate (sparse):")
  print(as.numeric(mrlasso.sparse)[1:25])
  print("l2 error:")
  print(l2.mrlasso.sparse)
  print("test data error:")
  print(error.mrlasso.sparse)
  
  ret.obj = list(m = m, n = n, d = d, noise = noise, s =  s,
                 l2.global = l2.global,
                 outlier.dist = outlier.dist,
                 error.global = error.global,
                 l2.adele.dense = l2.adele.dense,
                 l2.adele.sparse = l2.adele.sparse,
                 error.adele.dense = error.adele.dense,
                 error.adele.sparse = error.adele.sparse,
                 eta = eta,
                 l2.mrlasso.dense = l2.mrlasso.dense,
                 l2.mrlasso.sparse = l2.mrlasso.sparse,
                 error.mrlasso.dense = error.mrlasso.dense,
                 error.mrlasso.sparse = error.mrlasso.sparse)
  
  return(ret.obj)
}



# instance <- function(n = 100, m = 5, d = 2000, s = 5, 
#                      s.extra = 15, noise = 0.25, iter = 50)
# {
#   e = instance.i(n = n, m = m, d = d, s = s, noise = noise)
#   l2 = matrix(0, iter, length(e))
#   l2 = as.data.frame(l2)
#   colnames(l2) <- names(e)
#   l2[1, ] = unlist(e)
#   
#   for(i in 2:iter)
#   {
#     e = instance.i(n = n, m = m, d = d, s = s,
#                    noise = noise, s.extra = s.extra)
#     l2[i, ] = unlist(e)
#   }
#   return(l2)
# }


v.n = c(100 * 2^(0:4), rep(200, 5), rep(200, 5), rep(200, 5))
v.m = c(rep(5, 5), 2*2^(1:5), rep(5, 5), rep(5, 5))
v.s = c(rep(5, 5), rep(5, 5), 4 * 2^(0:4), rep(5, 5))
v.od = c(rep(25, 5), 10*2^(1:5), rep(25, 5), 25 * 2^(seq(-2, 2, by = 1)))

pars = data.frame(m = v.m, n = v.n, s = v.s, od = v.od)

n = 100; m = 8; d = 2000; s = 4; s.extra = 15;
noise = 0.25; snr = 4;
outlier.dist = 10; seed = 0;



args = commandArgs(trailingOnly=TRUE)
i = as.integer(args[1])


index = i %/% 50 + 1
m = as.integer(pars[index, 1])
n = as.integer(pars[index, 2])
s = as.integer(pars[index, 3])
od = as.numeric(pars[index, 4])


if(m == 5)
{
  s.extra = 20
} else {
  s.extra = 20
}
  
l2.list = instance.i(n = n, m = m, noise = 0.1, s = s, snr = 5,
                     s.extra = s.extra, seed = i+55999,
                     outlier.dist = od, d= 2000) 

write(rjson::toJSON(l2.list), file = paste('summaries_part/l2_', i, '.json', sep = ''))


