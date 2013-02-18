library(mvtnorm)
source("kernel.R")

N = 100
x = runif(N, 0, 1)
kappa = c(1,0.2,0)

Sigma = expSigma(kappa, x,x)
f = drop(rmvnorm(1, sigma=Sigma, method='svd'))
eps = rnorm(N,0,0.25)
y = f+eps

plot(x,y, pch=19)

marglike = function(kappa, y,x, sigma2)
{
	n = length(y)
	V = expSigma(kappa, x,x) + diag(sigma2, n)
	loglike = dmvnorm(y, sigma=V)
	return(loglike)
}

tau1grid = seq(0.01,5,length=50)
hgrid = seq(0.05,0.5,length=50)

kapgrid = as.matrix(expand.grid(tau1grid, hgrid))
kapgrid = cbind(kapgrid,0)
loglike = rep(0, nrow(kapgrid))
for(i in 1:nrow(kapgrid))
{
	loglike[i] = marglike(kapgrid[i,], y, x, 0.25^2)
}


loglike = matrix(loglike, nrow=length(tau1grid))
contour(tau1grid, hgrid, loglike)
points(kappa[1:2], pch=19, col='blue', cex=2)