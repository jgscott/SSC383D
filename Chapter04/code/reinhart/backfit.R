# Exercises 4
# Backfit an additive model to data.
# Alex Reinhart

source("poly-kernel-smoother.R")

# Return the kth partial residuals
# y: observations
# X: a matrix of predictor variables, with each row corresponding to one
#    observation in y and each column a new predictor
# k: column of X to ignore when computing residuals
# fs: prediction functions, corresponding to columns in X
# alpha: the fs are presumed to be fitted with the data mean-centered, so this
#        is the mean
kthPartials <- function(y, X, k, fs, alpha) {
  y - predictWith(as.matrix(X[,-k]), fs[-k], alpha)
}

# Given a matrix of data and a set of predictor functions, return estimated ys.
# X: a matrix of predictor variables, with each row corresponding to one
#    observation in y and each column a new predictor
# fs: a vector of prediction functions, corresponding to each column in X
# alpha: the fs are presumed to be fitted with the data mean-centered, so this
#        is the mean
predictWith <- function(X, fs, alpha) {
  y <- vector(mode="numeric", nrow(X))
  y[] <- alpha
  
  for (i in 1:ncol(X)) {
    if (!is.null(fs[[i]])) {
      y <- y + fs[[i]](X[,i])$val
    }
  }
  y
}

# Fit an additive model with backfitting. Uses cubic spline smoothing, since
# reuse of my kernel smoothing code caused inscrutable errors.
# y: a vector of responses
# X: a matrix of predictor variables, with each row corresponding to one
#    observation in y and each column a new predictor
# hs: bandwidths of kernel smoother for each variable
# tol: convergence parameter. Backfitting will continue until the fractional
#      change in the sum of squared errors is less than tol
backfit <- function(y, X, hs, tol) {
  alpha <- mean(y)
  
  fs <- vector("list", ncol(X))
  prevErr <- 1
  err <- 2
  while (abs(prevErr - err) / prevErr > tol) {
    for (j in 1:ncol(X)) {
      fs[[j]] <- fitPolySmoother(X[,j], kthPartials(y, X, j, fs, alpha),
                                 hs[j], gaussKernel)
    }
    prevErr <- err
    err <- sum((y - predictWith(X, fs, alpha))^2)
    print(err)
  }
  fs
}

air <- read.csv("../../examples/air.csv", header=TRUE)
X <- cbind(air$Wind, air$Temp, air$Solar.R)

# Fit ozone concentration vs. wind, temperature, and solar radiation
fs <- backfit(air$Ozone, X, c(0.5, 0.5, 10), 0.0001)

plot(air$Wind, air$Ozone)
xs = seq(from=min(air$Wind), to=max(air$Wind))
lines(xs, fs[[1]](xs)$val + mean(air$Ozone))

plot(air$Temp, air$Ozone)
xs = seq(from=min(air$Temp), to=max(air$Temp))
lines(xs, fs[[2]](xs)$val + mean(air$Ozone))

plot(air$Solar.R, air$Ozone)
xs = seq(from=min(air$Solar.R), to=max(air$Solar.R))
lines(xs, fs[[3]](xs)$val + mean(air$Ozone))
