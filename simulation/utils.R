l2.norm <- function(x)
{
  ind = which(x != 0)
  return(sqrt(sum(x[ind] ^2 )))
}


l.infty.norm <- function(x)
{
  return(max(abs(x)))
}


soft.th <- function(beta, lambda = 0.1)
{
  non_zeros = which(beta != 0)
  for(i in non_zeros)
  {
    z = beta[i]
    beta[i] = sign(z) * max(c(abs(z) - lambda, 0))
  }
  return(beta)
}



eval.central.beta <- function(data.test, beta.central)
{
  s = 0
  count = 0
  x = data.test$x
  y = data.test$y
  
  for(i in 1:length(x))
  {
    x.i = x[[i]]
    y.i = y[[i]]
    count = count + length(y.i)
    error.i = y.i - x.i %*% beta.central[-1] - beta.central[1]
    s = s + sum(error.i^2)
  }
  return(s/count)
}
