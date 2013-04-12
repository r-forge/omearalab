if(!require(gmp)) {
  install.packages("gmp",repos="http://cran.us.r-project.org", lib=tempdir())
  .libPaths(c(.libPaths(), tempdir()))
}
source("V7_UtilityFns.R")
source("V7_NewSimulator.R")
source("V7_StochasticSSASims_Functions.R")
load("Rates.Rsave")