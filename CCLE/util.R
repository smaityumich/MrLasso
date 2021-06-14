soft.th <- function(x, lambda = 0.1)
{
  r = rep(0, length(x))
  non_zeros = which(x != 0)
  for(i in non_zeros)
  {
    z = x[i]
    r[i] = sign(z) * max(c(abs(z) - lambda, 0))
  }
  return(r)
}