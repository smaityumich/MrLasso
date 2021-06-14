library(Matrix)

get.beta <- function(m, d = 10000, s = 1)
{
  beta = Matrix(0, d, m, sparse = T)
  beta[1:2,] = s
 
  for(k in 1:m)
    beta[k+2, k] = s * 8
  return(beta)
}


get.beta0 <- function(d = 10000, s = 1)
{
  beta = Matrix(0, d+1, 1, sparse = T)
  beta[2:3] = s
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
    ind = c(1, 2, k+2)
    y[[k]] = x[[k]][, ind] %*% beta[ind, k] + rnorm(n) * noise
  }
  return(y)
}


get.data <- function(m, d, n, s = 1, noise = 0.05)
{
  p = list(m = m, n = n, d = d, s = s, noise = noise)
  beta = get.beta(m, d, s = s)
  beta.mrlasso = get.beta0(d, s= s)
  beta.adele = beta.mrlasso
  beta.adele[(3+1):(3+m)] = s * 8/m
  
  
  x = get.x(n, m, d)
  y = get.y(beta, x, noise = noise)
  return(list(x = x, y = y, beta.adele = beta.adele, beta.mrlasso = beta.mrlasso, pars = p, beta = beta))
}


