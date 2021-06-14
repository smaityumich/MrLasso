library(dplyr)
library(glmnet)
library(EnvStats)
library(gridExtra)
setwd('D:/MRLASSO/CCLE_modified/')
source('utils_plot.R')


folder = 'ccle_cv/'

drug.data = read.csv('CCLE_NP24.2009_Drug_data_2015.02.24.csv')
drugs = unique(drug.data$Compound)

for(drug in drugs)
{
  file = paste(folder, 'summary_', drug, '.RData',  sep = '')
  if(file.exists(file))
  {
    load(file = paste(folder, 'summary_', drug, '.RData',  sep = ''))
    
    
    p.majority = cat.boxplot(cv.global,   cv.adele.best, cv.mrlasso.best, total,
                             cat = 'majority', legend.position = 'none', 
                             x.text.angle = 40, prune = 6, 
                             ylim = c(0, 1))
    
    p.minority = cat.boxplot(cv.global,   cv.adele.best, cv.mrlasso.best, total,
                             cat = 'minority', legend.position = 'bottom', 
                             x.text.angle = 40, prune = 6, 
                             ylim = c(0, 1))
    
    
    coef.plot = coeff.plot(betas.local)
    
    
    organs.cv = grid.arrange(p.majority, p.minority, ncol = 1, heights = c(1, 1.2))
    
    pdf(paste(folder , drug, '.pdf', sep = ''), height = 6, width = 9)
    grid.arrange(organs.cv, coef.plot, nrow = 1, widths = c(1.8, 0.8),
                 top = textGrob(paste('Drug:', drug))
    )
    dev.off()
  }
  
  
}
