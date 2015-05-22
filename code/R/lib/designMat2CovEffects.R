designMat2CovEffects <-function(COVARIATES_MAP, sigCovars.Effects){
  sapply(COVARIATES_MAP,function(x,y = sigCovars.Effects){return(mean(y[x]))})
}