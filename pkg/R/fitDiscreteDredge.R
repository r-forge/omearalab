library(parallel)
library(geiger)

fitDiscreteTransformFirst(x, phy, data) {
  return(fitDiscrete(phy=phy, data=data, treeTransform=x))
}

fitDiscreteDredge<-function(phy, data, mc.cores=detectCores()) {
  transforms<-c("none", "lambda", "kappa", "delta", "linearChange", "exponentialChange", "twoRate")
  results<-mclapply(transforms, fitDiscreteTransformFirst, phy=phy, data=data, mc.cores=mc.cores)
  return(results)
}
