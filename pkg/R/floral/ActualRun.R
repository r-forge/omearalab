if(!require(gmp)) {
  install.packages("gmp",repos="http://cran.us.r-project.org", lib=tempdir())
  .libPaths(c(.libPaths(), tempdir()))
}
source("V6_UtilityFns.R")
source("V6_NewSimulator.R")
source("V6_StochasticSSASims_Functions.R")
load("Rates.Rsave")

doParallelSSA(tf=136, x0=x0, q.means=q.means, lambda.means=lambda.means, mu.means=mu.means, maxWallTime=Inf, file.string="FullHistory", rescale.species=250000, yule.scale=0, full.history=FALSE,print.freq=10000)
