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

