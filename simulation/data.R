library(Matrix)
source('mrlasso.R')


beta.values.1 <- function(mu, delta, eps = 0, sd = 0.1)
{
  z = rbinom(1, 1, 1-eps)
  if(z)
    return(rnorm(1, mean = mu, sd = sd))
  else
  {
    return(runif(1, mu + delta, mu + delta + 1))
  }
}


beta.values.2 <- function(delta)
{
  r = 0
  if(rbinom(1, 1, 0.5))
    r = runif(1, delta, delta + 1)
  return(r)
}







get.beta <- function(m, d = 10000, s = 5, s.extra = 15, snr = 1, outlier.dist = 10, seed = 0)
{
  set.seed(seed)
  BETA = matrix(0, d, m)
  for(i in 1:s)
  {
    for(j in 1:m)
    {
      BETA[i, j] = beta.values.1(mu = snr, delta = outlier.dist, eps = 0)
    }
  }
  for(i in (s+1):(s + s.extra))
  {
    j = sample.int(m, size = 1)
    BETA[i, j] = beta.values.2(delta = snr + outlier.dist)
  }
  BETA = Matrix(BETA, sparse = T)

  
  return(BETA)
}


get.beta.test <- function(m, d = 10000, s = 5, snr = 1, outlier.dist = 10, seed = 0)
{
  set.seed(seed)
  BETA = matrix(0, d, m)
  for(i in 1:s)
  {
    for(j in 1:m)
    {
      BETA[i, j] = beta.values.1(mu = snr, delta = outlier.dist, eps = 0, sd = 0)
    }
  }
  
  BETA = Matrix(BETA, sparse = T)
  
  
  return(BETA)
}


get.beta0.mrlasso <- function(BETA, s = 5, s.extra = 15)
{
  d = BETA@Dim[1]
  beta = Matrix(0, d+1, 1, sparse = T)
  beta[2:(s + s.extra + 1)] = apply(BETA[1:(s + s.extra), ], 1, function(x){return(mr.lasso(x, eta = 3))})
  return(beta)
}

get.beta0.adele <- function(BETA)
{
  d = BETA@Dim[1]
  beta = Matrix(0, d+1, 1, sparse = T)
  beta[2:(d+1)] = apply(BETA, 1, mean)
  return(beta)
}



get.x <- function(n, m, d = 10000)
{
  x = list()
  for(k in 1:m)
  {
    x.k = matrix(0, n, d)
    for(i in 1:n)
    {
      x.k[i, ] = rnorm(d)
    }
    x[[k]] = x.k
  }
  return(x)
}



get.y <- function(beta, x, noise = 0.05)
{
  m = ncol(beta)
  d = nrow(beta)
  n = nrow(x[[1]])
  
  y = list()
  for(k in 1:m)
  {
    y[[k]] = x[[k]] %*% beta[, k] + rnorm(n) * noise
  }
  return(y)
}


get.data <- function(m = 5, d = 100, n = 20, snr = 1, noise = 0.05, s = 5, s.extra = 15, outlier.dist = 10, seed = seed)
{
  p = list(m = m, n = n, d = d, s = s, s.extra = s.extra,
           noise = noise, snr = snr, outlier.dist = outlier.dist)
  beta = get.beta(m, d = d, s = s, s.extra = s.extra,
                  snr = snr, outlier.dist = outlier.dist, seed = seed)
  beta.mrlasso = get.beta0.mrlasso(beta, s = s, s.extra = s.extra)
  beta.adele = get.beta0.adele(beta)
  
  
  x = get.x(n, m, d)
  y = get.y(beta, x, noise = noise)
  return(list(x = x, y = y, beta.adele = beta.adele, beta.mrlasso = beta.mrlasso, pars = p, beta = beta))
}


get.data.test <- function(m = 5, d = 100, n = 20, snr = 1, noise = 0.05, s = 5, s.extra = 15, outlier.dist = 10, seed = seed)
{
  p = list(m = m, n = n, d = d, s = s, s.extra = s.extra,
           noise = noise, snr = snr, outlier.dist = outlier.dist)
  beta = get.beta.test(m, d = d, s = s,
                  snr = snr, outlier.dist = outlier.dist, seed = seed)
  beta.mrlasso = get.beta0.mrlasso(beta, s = s, s.extra = s.extra)
  beta.adele = get.beta0.adele(beta)
  
  
  x = get.x(n, m, d)
  y = get.y(beta, x, noise = noise)
  return(list(x = x, y = y, beta.adele = beta.adele, beta.mrlasso = beta.mrlasso, pars = p, beta = beta))
}

