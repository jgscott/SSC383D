library(Rcpp)
library(MASS)
sourceCpp('normix.cpp')

weights=c(0.5,rep(0.1,5))
mu=c(0, (0:4)/2 - 1)
tau2 = c(1, rep(1/100, 5))

# Density function
curve(dnormix(x, weights, mu, tau2), from=-3, to=3, n=1001)

# Draw a sample
N = 5000
X = rnormix(N, weights, mu, tau2)
hist(X, 100)

par(mfrow=c(1,2))
d1 = density(X, bw=1, kernel='epanechnikov')
truehist(X, 100, main='1', ylim=c(0,0.75))
lines(d1, lwd=3)

d2 = density(X, bw=.01, kernel='epanechnikov')
truehist(X, 100, main='.005', ylim=c(0,0.75))
lines(d2, lwd=3)


par(mfrow=c(1,2))
d3 = density(X, bw=.1, kernel='epanechnikov')
truehist(X, 100, main='0.1', ylim=c(0,0.75))
lines(d3)

d4 = density(X, bw='bcv', kernel='epanechnikov')
truehist(X, 100, main='Cross Validation', ylim=c(0,0.75))
lines(d4)

