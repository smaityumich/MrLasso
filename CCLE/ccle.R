library(dplyr)
library(glmnet)
library(EnvStats)
library(gridExtra)


## Seed for reproducibility
set.seed(0)

setwd('D:/MRLASSO/CCLE_modified/')
source('utils_ccle.R')
source('utils_lasso.R')
source('utils_adele_ccle.R')
source('utils_mrlasso_ccle.R')
source('utils_final.R')
source('utils_plot.R')

folder = 'ccle_cv/'

if(!dir.exists(folder))
  dir.create(folder)


## Data download portal: https://depmap.org/portal/download/
## DepMap Public 21Q2 data

cell.line = read.csv('CCLE_expression.csv')
drug.data = read.csv('CCLE_NP24.2009_Drug_data_2015.02.24.csv')
cell.line.sep = strsplit(as.character(drug.data$CCLE.Cell.Line.Name), split = '_')

organ = rep('', length(cell.line.sep))
for(i in 1:length(cell.line.sep))
{
  organ.split = cell.line.sep[[i]][-1]
  
  if(length(organ.split) == 0)
    organ[i] = NA
  else
    organ[i] = do.call(paste, c(as.list(organ.split), sep = "_"))
}
drug.data$Organ = organ
drug.data = drug.data[which(drug.data$FitType == 'Sigmoid'), ]
rm(cell.line.sep)

secondary.dose = read.csv('secondary-screen-dose-response-curve-parameters.csv')
id.pair = subset(secondary.dose, select = c(depmap_id, ccle_name))
rm(secondary.dose)
id.pair = id.pair %>% distinct()



drugs = unique(drug.data$Compound)







for(i.drug in 1:length(drugs))
{
  drug = drugs[i.drug]
  ct = 10
  total = f.drug.data(drug, drug.data, cell.line, id.pair)
 
    
    
    
    
    
    organ.count = table(total$Organ)
    majority.organ = organ.count[which(organ.count>ct)]
    majority.organs = names(majority.organ)
    
    if(length(majority.organs) > 2)
    {
      minority.organ = organ.count[which(organ.count<=ct)]
      minority.organs = names(minority.organ)
      
      test.index.majority = test.index.organs(total, majority.organs, save.index = F, folds = 10)
      index.minority = get.minority.indices(total, minority.organs)
      
      majority.total.index = which(total$Organ %in% majority.organs)
      
      
      
      # Global lasso
      cv.global = global.lasso(total, majority.total.index, test.index.majority, index.minority)
      cv.mean = global.mean(total, majority.total.index, test.index.majority, index.minority)
      
      beta.cv = local.lasso(total, test.index.majority)
      
      
      
      logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
      lambdas = logspace(-1.5, 0.5, 10)
      etas = logspace(-1.5, 0.5, 10)
      
      
      adele.cvs.lambdas = rep(0, length(lambdas))
      for(i.lambda in 1:length(lambdas))
      {
        lambda = lambdas[i.lambda]
        cv.adele = mse.adele(beta.cv, total, test.index.majority, index.minority, lambda = lambda)
        adele.cvs.lambdas[i.lambda] = summary.mse(cv.adele, total)
      }
      
      
      lambda.adele.best = min(adele.cvs.lambdas)
      
      # ADELE with best CV-mse
      cv.adele.best = mse.adele(beta.cv, total, test.index.majority, index.minority, lambda = lambda.adele.best)
      
      
      
      beta.mrlasso.debiased.eta.all = mr.lasso.all(beta.cv = beta.cv, etas = etas)
      
      mrlasso.cv.lambda.eta = matrix(0, length(etas), length(lambdas))
      
      for(i.eta in 1:length(etas))
      {
        obj = beta.mrlasso.debiased.eta.all$betas[[i.eta]]
        for(i.lambda in 1:length(lambdas))
        {
          lambda = lambdas[i.lambda]
          cv.mrlasso = cv.integrated.mrlasso.eta.lambda(obj, total, test.index.majority,
                                                        index.minority, lambda)
          mrlasso.cv.lambda.eta[i.eta, i.lambda] = summary.mse(cv.mrlasso, total)
        }
        
      }
      
      i.best = which(mrlasso.cv.lambda.eta==min(mrlasso.cv.lambda.eta),arr.ind=TRUE)
      eta.mrlasso.best = etas[i.best[1]]
      lambda.best.mrlasso = lambdas[i.best[2]]
      
      
      
      
      obj.best = beta.mrlasso.debiased.eta.all$betas[[i.best[1]]]
      # MRLASSO with best CV-mse
      cv.mrlasso.best = cv.integrated.mrlasso.eta.lambda(obj.best, total, test.index.majority,
                                                         index.minority, lambda.best.mrlasso)
      
      betas.local = local.lasso.final(total, majority.organs, lambda.adele.best, lambda.best.mrlasso, eta.mrlasso.best)
      
      
      rm(list = c('beta.mrlasso.debiased.eta.all', 'beta.cv', 'cv.adele', 'cv.mrlasso', 'obj', 'obj.best'))
      
      save.list = c('betas.local', 'total', 'cv.global', 'cv.mean', 'cv.adele.best',
                    'cv.mrlasso.best', 'index.minority', 'test.index.majority',
                   'lambda.adele.best',
                    'lambda.best.mrlasso', 'eta.mrlasso.best', 'etas', 'lambdas',
                    'majority.organs', 'minority.organs')
      
      save(list = save.list, file = paste(folder, 'summary_', drug, '.RData',  sep = ''))
      
      
    }
  
}



