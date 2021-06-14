library(dplyr)
library(glmnet)
library(ggplot2)
library(stringr)
library(tidyverse)
library(latex2exp)
library(scales)
library(grid)



overall.cv <- function(cv, p.organs)
{
  w.mean <- function(x)
  {
    return(sum(x*p.organs))
  }
  overall = apply(cv, 2, w.mean)
  return(rbind(cv, overall))
}

modify.str <- function(x)
{
  for(i in 1:length(x))
  {
    organ = x[i]
    if(organ == 'CENTRAL_NERVOUS_SYSTEM')
      organ.modified = 'CNS'
    else if(organ == 'UPPER_AERODIGESTIVE_TRACT')
      organ.modified = 'UAT'
    else
      organ.modified = organ
    organ.modified = str_replace_all(organ.modified, '_', ' ')
    organ.modified = tolower(organ.modified)
    x[i] = organ.modified
  }
  return(x)
}



modify.mat <- function(cv, name = 'global')
{
  cv = data.frame(cv)
  colnames(cv) = 1:ncol(cv)
  df = cv %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
  colnames(df) <- c('organ', 'i.cv', 'mse')
  df[['method']] = name
  return(df)
}

cat.boxplot <- function(cv.global, cv.adele.best, cv.mrlasso.best, total,
                        cat = 'majoirty', legend.position = 'none', 
                        x.text.angle = 40, prune = 4, 
                        ylim = c(0, 1.5))
{
  
  organs = row.names(cv.global[[cat]])
  n.organs = rep(0, length(organs))
  for(i.organ in 1:length(organs))
    n.organs[i.organ] = length(which(total$Organ == organs[i.organ]))
  
  p.organs = n.organs/sum(n.organs)
  
  global.majority = overall.cv(cv.global[[cat]], p.organs)
  mrlasso.majority = overall.cv(cv.mrlasso.best[[cat]], p.organs)
  adele.majoirty = overall.cv(cv.adele.best[[cat]], p.organs)
  
  df = rbind(
    modify.mat(global.majority, 'global'),
    modify.mat(mrlasso.majority, 'MR.Lasso'),
    modify.mat(adele.majoirty, 'ADELE')
  )
  
  
  i.pruned.organs = which(n.organs > prune)
  pruned.organs = c(organs[i.pruned.organs], 'overall')
  df.pruned = df[df$organ %in% pruned.organs,]
  
  df.pruned$organ = modify.str(df.pruned$organ)
  
  orgs <- unique(df.pruned$organ)
  orgs <- orgs[which(orgs != 'overall')]
  orgs = c(orgs, 'overall')
  
  df.pruned$organ = factor(df.pruned$organ, labels = orgs)
  
  p <-ggplot(df.pruned, aes(x=organ, y=mse, fill=method)) +
    geom_boxplot(width = 0.6, lwd = 0.1, outlier.shape = NA) + 
    theme_bw() + 
    scale_fill_discrete(labels = c("ADELE", "Mr.Lasso", 'global.Lasso'), name = '')  + 
    theme(axis.text.x = element_text(angle=x.text.angle, hjust=1), legend.position = legend.position) +
    xlab(element_blank()) + 
   # ylim(c(0, 10)) + 
    ylab('MSE') + 
    scale_y_continuous(trans = 'log2', breaks = 2^(seq(-6, 4, by = 2)),
                       labels = trans_format("log2", math_format(2^.x)))
  
  
  return(p)
}




coeff.plot <- function(betas.local)
{
  i.nzeros = c()
  #for(organ in betas.local$beta)
  #  i.nzeros = c(i.nzeros, (which(organ[['beta']][-1] != 0) + 1))
  i.nzeros = c(i.nzeros, (which(betas.local$adele[-1] != 0) + 1))
  i.nzeros = c(i.nzeros, (which(betas.local$mrlasso[-1] != 0) + 1))
  #i.nzeros = c(i.nzeros, (which(betas.local$global.lasso$beta[-1] != 0) + 1))
  i.nzeros = unique(i.nzeros)
  i.nzeros = sort(i.nzeros)
  
  name.coeff = row.names(betas.local$adele)
  
  if(length(i.nzeros) > 0)
  {
    df.coeffs = data.frame(ADELE = betas.local$adele[i.nzeros],Mr.Lasso = betas.local$mrlasso[i.nzeros])#,
                           #global.Lasso = betas.local$global.lasso$beta[i.nzeros])
    
    
    organs = names(betas.local$beta)
    for(organ in organs)
    {
      if(organ == 'CENTRAL_NERVOUS_SYSTEM')
        organ.modified = 'CNS'
      else if(organ == 'UPPER_AERODIGESTIVE_TRACT')
        organ.modified = 'UAT'
      else
        organ.modified = organ
      organ.modified = str_replace_all(organ.modified, '_', ' ')
      organ.modified = tolower(organ.modified)
      df.coeffs[[organ.modified]] = betas.local$beta[[organ]][['beta']][i.nzeros]
    }
    rownames(df.coeffs) = name.coeff[i.nzeros]
    
    dt <- df.coeffs %>%
      rownames_to_column() %>%
      gather(colname, value, -rowname)
    head(dt)
    dt[['size']] = abs(dt$value)/max(abs(dt$value))
    
    color = rep('blue', nrow(dt))
    for(i in 1:nrow(dt))
    {
      if(dt$value[i] < 0)
        color[i] = 'red'
    }
    dt[['color']] = color
    #dt$rowname = str_remove(dt$rowname, 'ENSG00000')
    
    dt$rowname = str_replace(str_split(dt$rowname, '\\.\\.', simplify = T)[, 1], '\\.', '-')
    
    
    n.adele = sum(df.coeffs$ADELE != 0)
    n.mr.lasso = sum(df.coeffs$Mr.Lasso != 0)
    
    if(T)#(n.mr.lasso > 0)
    {
      level_order <- c(colnames(df.coeffs)[-(1:2)], 'ADELE', 'Mr.Lasso')#, 'global.Lasso')
      
      #name.order <- c(colnames(df.coeffs)[-c(1,2)],  paste('ADELE\n(s=', n.adele, ')', sep = ''),
      #                 paste('Mr.Lasso\n(s=', n.mr.lasso, ')', sep = ''))
      name.order <- c(colnames(df.coeffs)[-(1:2)],  'ADELE', 'Mr.Lasso')#, 'global.Lasso')
      dt_part = dt[dt$value != 0, ]
      coef.plot = ggplot(dt_part, aes(x = factor(colname, levels = level_order), y = rowname, fill = color, color = color, size = size))+
        geom_point() + 
        theme_bw() +
        scale_color_manual(breaks = c("red", "blue"),
                           values=c("#F8766D", "#00BFC4")) + 
        theme(axis.text.y = element_text(size = 7)) +
        theme(axis.text.x = element_text(angle=45, hjust=1, size = 10), legend.position = 'none') + 
        xlab(element_blank()) + 
        ylab(element_blank()) + 
        scale_x_discrete(labels = level_order, limits = c(colnames(df.coeffs)[-(1:2)],  'ADELE',  'Mr.Lasso')) + #,
                                                         #'global.Lasso')) +
        guides(size = F)
      return(coef.plot)
    }
    else
      return(NULL)
  }
  else
    return(NULL)
}