library(mvtnorm)

rtnorm = function(n, mu, sigma, lower=-Inf, upper=Inf)
# Simulated truncated normals using the inverse CDF
{
	u.l = pnorm(lower, mu, sigma)
	u.u = pnorm(upper, mu, sigma)
	u = runif(n, u.l, u.u)
	z = qnorm(u, mu, sigma)
	z;
}

polls = read.csv("polls.csv", header=TRUE)

polls = na.omit(polls)

N = nrow(polls)

# Extract the data and create the design matrices for fixed and random effects
Y = 0+{polls$bush==1}
X.r = model.matrix(bush ~ state-1, data=polls)
X.f = model.matrix(bush ~ edu + black + female, data=polls)[, 2:6]
P.r = ncol(X.r)
P.f = ncol(X.f)

# Priors 
# tau2.theta ~ IG(a/2,b/2)
a = 1
b = 1
tau2.beta = 1e6

# Initialize
theta = rep(0, ncol(X.r))
beta = rep(0, ncol(X.f))
tau2.theta = 1

# The bounds for z implied by Y
zl = rep(-Inf,N)
zu = rep(Inf, N)
zl[Y==1] = 0
zu[Y==0] = 0

burn = 500
NMC = 2000
betasave = matrix(0, nrow=NMC, ncol=P.f)
thetasave = matrix(0, nrow=NMC, ncol=P.r)
for(nmc in 1:(NMC+burn))
{
	if(nmc %% 50 == 0) cat("Iteration", nmc, "\n")
	# First sample the latent z's
	zhat = X.r %*% theta + X.f %*% beta
	z = rtnorm(N, zhat, 1, zl, zu)
	
	# Now the fixed effects
	zpartial = z - X.r %*% theta
	beta.Sig = solve(crossprod(X.f) + diag(1/tau2.beta, P.f))
	beta.mu = beta.Sig %*% crossprod(X.f, zpartial)
	beta = drop(rmvnorm(1, beta.mu, beta.Sig))
	
	# Now the random effects
	zpartial = z - X.f %*% beta
	theta.Sig = solve(crossprod(X.r) + diag(1/tau2.theta, P.r))
	theta.mu = theta.Sig %*% crossprod(X.r, zpartial)
	theta = drop(rmvnorm(1, theta.mu, theta.Sig))
	
	# Now the hyperparameters
	bnew = b + sum(theta^2)
	tau2.theta = 1/rgamma(1, a + P.r, rate=bnew)
	
	# Save the draws
	if(nmc > burn)
	{
		betasave[nmc-burn,] = beta
		thetasave[nmc-burn,] = theta
	}
}

colnames(thetasave) = levels(polls$state)
thetahat = colMeans(thetasave)
myorder = order(thetahat)

boxplot(thetasave[,myorder], las=2)
boxplot(pnorm(thetasave[,myorder]), las=2)

# Compute the raw fractions of Bush respondents in each state
library(mosaic)
statesize = summary(polls$state)
rawmean = mean(bush~state, data=polls)

# Compare these versus the posterior means
# Not exactly comparable b/c of covariates, but close
plot(rawmean, pnorm(thetahat))
abline(0,1)

plot(statesize, rawmean)
plot(statesize, thetahat)

# Plot the amount of shrinkage versus state-level sample size
plot(statesize, (pnorm(thetahat)-rawmean)/rawmean)

