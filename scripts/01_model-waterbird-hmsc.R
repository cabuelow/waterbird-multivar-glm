# joint species distribution model (JSDM)

# libaries

library(Hmsc)
library(tidyverse)
library(GGally)

# data

birds.env <- read.csv('data/waterbird-basin-threats.csv')
preds <- birds.env %>% dplyr::select(Consumptive_Water_Loss, # select predictors only
                                     Pesticide_Loading, 
                                     Phosphorus_Loading, Nitrogen_Loading,
                                     Aquaculture_Pressure, Mercury_Deposition, 
                                     ire_pc_sse, pre_mm_syr)

# check multi-collinearity in predictor variables

ggcorr(preds, label = TRUE)

# specify hmsc model

xycoords <- as.matrix(birds.env %>% dplyr::select(X,Y))
colnames(xycoords) <- c('x-coordinate', 'y-coordinate')
row.names(xycoords) <- 1:nrow(birds.env)
y <- as.matrix(birds.env %>% dplyr::select(Australasian.Darter:Australasian.Bittern))
x <- preds

# spatial random effect

studyDesign <- data.frame(sample = as.factor(1:nrow(birds.env)))
rL.spatial <- HmscRandomLevel(sData = xycoords)
rL.spatial <- setPriors(rL.spatial, nfMin=3, nfMax=5) # here set number of LVs (3-5)

m <- Hmsc(Y = y, XData = x, 
          XFormula = as.formula(paste('~', paste(colnames(x), collapse = '+'))),
          studyDesign = studyDesign, 
          ranLevels = list(sample = rL.spatial), 
          distr = 'probit')

# set mcmc sampling parameters 

nChains <- 4
test.run <- TRUE
if (test.run){
  thin = 1
  samples = 10
  transient = 5
  verbose = 0
} else {
  thin = 100
  samples = 1000
  transient = ceiling(0.5*samples*thin)
  verbose = 1 # set to 0 to suppress output reporting how MCMC sampling is proceeding
}

# run model

set.seed(123)

system.time( # takes ~4 hr to run
m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient, adaptNf = rep(ceiling(0.4*samples*thin),1),
                       nChains = nChains, nParallel = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
)

# save model

saveRDS(m, 'outputs/models/mod-spatialRF.rds')



