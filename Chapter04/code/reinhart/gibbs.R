# Exercises 4
# Use Gibbs sampling to fit a model to data.
# Alex Reinhart

library(mvtnorm)

# Perform Gibbs sampling to determine the fit matrix B in a Bayesian multiple-
# regression model, using the given hyperparameters for the priors.
# X: a design matrix
# y: response vector
# N: number of Gibbs iterations to run. Averages over results for last quarter.
# Returns B, a vector whose elements are regression coefficients corresponding
# to columns of X.
gibbsFit <- function(X, y, N, a, b, c, d) {
  # Setup
  n <- nrow(X)
  p <- ncol(X)
  Br <- matrix(nrow=p, ncol=N)

  # Create initial draws from prior.
  sigma2 <- 1/rgamma(1, a/2, rate=b/2)
  tau2 <- 1/rgamma(1, c/2, rate=d/2)
  
  for (i in 1:N) {
    mn <- solve(((1/sigma2) * t(X) %*% X + 1/tau2)) %*% ((1/sigma2) * t(X) %*% y)
    sd <- solve(1/sigma2 * t(X) %*% X + 1/tau2)
    B <- t(as.matrix(rmvnorm(1, mean=mn, sigma=sd)))
    Br[,i] <- B
    
    sigma2 <- 1/rgamma(1, (n+a)/2, rate=(b + t(y - X %*% B) %*% (y - X %*% B)) / 2)
    tau2 <- 1/rgamma(1, (p+c)/2, rate=(t(B) %*% B + d)/2)
  }
  rowMeans(Br[,floor(N/4):N])
}
