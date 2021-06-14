library(glmnet)


f.drug.data <- function(drug, drug.data, cell.line, id.pair)
{
  drug.data.1 = drug.data[which(drug.data$Compound == drug),]
  drug.data.1 = subset(drug.data.1, select = c(CCLE.Cell.Line.Name, ActArea, Organ, Compound))
  colnames(drug.data.1)[1] = 'ccle_name'
  
  
  
  cell.line.1 = cell.line[cell.line$X %in% id.pair$depmap_id, ]
  colnames(cell.line.1)[1] = 'depmap_id'
  
  
  total = merge(id.pair, drug.data.1, by = 'ccle_name')
  total = merge(total, cell.line.1, by = 'depmap_id')
  rm(list = c('cell.line.1', 'drug.data.1'))
  total = total[complete.cases(total), ]
  total = subset(total, select = -c(depmap_id, ccle_name, Compound))
  return(total)
}





test.index <- function(index, folds = 3, seed = 0)
{
  set.seed(seed)
  #ind.randomized = sample(index, length(index), replace = F)
  s.fold = floor(length(index)/5) # size of each fold
  r.index = matrix(0, folds, s.fold)
  for(i in 1:folds)
    r.index[i, ] = sample(index, s.fold, replace = F) #ind.randomized[((i-1) * s.fold + 1):(i*s.fold)]
  return(r.index)
}



test.index.organs <- function(total, organs, folds = 3, seed = 0,
                                       filename = '.RData', load.index = F, save.index = T)
{
  if(load.index)
  {
    load(file = filename)
  }
  else
  {
    test.index.majority = list()
    for(organ in organs)
    {
      index = which(total$Organ == organ)
      test.index.majority[[organ]] = test.index(index, folds = folds, seed = seed)
    }
    if(save.index)
      save(list = c('test.index.majority'), file = filename)
  }
  return(test.index.majority)
}


get.minority.indices <- function(total, organs)
{
  r = list()
  for(organ in organs)
  {
    r[[organ]] = which(total$Organ == organ)
  }
  return(r)
}




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



#lp= local.lasso(total, test.index.majority)
#obj = lp[[1]]





summary.mse <- function(return.mse, total)
{
  
  mse = return.mse[['majority']]
  organs = rownames(mse)
  
  n.organs = rep(0, length(organs))
  for(i.organ in 1:length(organs))
    n.organs[i.organ] = length(which(total$Organ == organs[i.organ]))
  
  p.organs = n.organs/sum(n.organs)
  mean.tr <- function(x)
  {
    i = which.max(x)
    return(mean(x[-i]))
  }
  mse.organs = apply(mse, 1, mean)
  return(sum(p.organs * mse.organs))
}
