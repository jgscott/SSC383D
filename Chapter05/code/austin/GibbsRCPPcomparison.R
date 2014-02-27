library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)

setwd("C:/Users/Austin/Dropbox/Stat Mod 2/Exercise 4 R scripts")

sourceCpp('rcppgibbsamp.cpp')

perform <- benchmark(gibbs(10000,ydc,Xd,1,1,1,1,1,1),gibbsC(10000,ydc,Xd,1,1,1,1,1,1),replications=1)