# Sparse matrices in the Eigen library
library(Matrix)
library(RcppEigen)
sourceCpp("sgdlogitsparse.cpp")

# Simulate a biggish data set
N = 1e7
P = 500
rx = c(1:N,sample(1:N, 9*N, replace=TRUE))
cx = sample(1:P, 10*N, replace=TRUE)
xx = rnorm(length(rx))

# store with cases in column-major format
X = sparseMatrix(i = cx, j = rx, x=xx)
beta = rnorm(P)
psi = drop(crossprod(X, beta))
w = 1/{1+exp(-psi)}
y = rbinom(N, 1, w)
n = as.integer(rep(1,N))

system.time(sgd1 <- sparsesgd_logit(X, y, n, pars=c(10,.51), npass=5))
betahat = sgd1$beta

plot(beta, betahat); abline(0,1)

# try if you dare
system.time(glm1 <- glm(y ~ t(as.matrix(X))-1, family=binomial))

plot(beta, coef(glm1))

plot(coef(glm1), betahat)