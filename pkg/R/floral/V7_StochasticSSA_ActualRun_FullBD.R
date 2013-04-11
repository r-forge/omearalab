if(!require(gmp)) {
  install.packages("gmp",repos="http://cran.us.r-project.org", lib=tempdir())
  .libPaths(c(.libPaths(), tempdir()))
}
source("V7_UtilityFns.R")
source("V7_NewSimulator.R")
source("V7_StochasticSSASims_Functions.R")
load("Rates.Rsave")

#tf is based on age of angiosperms (136) plus 20 MY in future
#t.rescale is based on age of angiosperms and current diversity
doParallelSSA(tf=156, x0=x0, q.means=q.means, lambda.means=lambda.means, mu.means=mu.means, maxWallTime=Inf, file.string="FullHistory_bd", rescale.species=250000, yule.scale=0, full.history=FALSE,print.freq=10000, t.rescale=136)
