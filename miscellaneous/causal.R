tstat = function(mylm, nullval=0) {
	return({coef(mylm) - nullval}/sqrt(diag(vcov(mylm))))
}

# Set the beta (the true X on Y effect) and N (sample size)
BetaYonX = 1
N = 50

# Simulate the data from a structural model
Z = rnorm(N, 0, 1)
X = Z + rnorm(N, 0, 1)
A = 3 + 2*Z + rnorm(N, 0, 1)
B = -3 + 2*Z + rnorm(N, 0, 1)
Y = BetaYonX*X + A + B + rnorm(N, 0, 1)
D = 2 + X + Y + rnorm(N, 0, 1) 

# Ignoring confounding
lm0 = lm(Y~X)
tstat(lm0, 1)

# Fit the model using various "back-door" strategies
lm1 = lm(Y ~ X + Z)
lm2 = lm(Y ~ X + A + B)

# Now conditioning on a descendant
lm3 = lm(Y ~ X + Z + D)




####
# Front door adjustment

N = 1000

# Our structural model
U = rnorm(N, 0, 1)  # unobserved common cause of X and Y
X = U + rnorm(N, 0, 1)  # covariate we care about
M = X + rnorm(N, 0, 1)  # mediating variable ("front door path")
Y = M + U + rnorm(N, 0, 1) # outcome

# An inconsistent estimate
lm1 = lm(Y~X)
summary(lm1)

# First estimate the effect of X on M
lm2 = lm(M~X)
summary(lm2)

# Then the effect of M on Y
# Need to condition on X to block the correlating back-door path via U
# Notice the X coefficient in this model is the wrong causal estimate
lm3 = lm(Y~M+X)
summary(lm3)

# Compose to get a consistent estimate of causal effect of X on Y
# without knowledge of the unobserved common cause U
coef(lm2)[2] * coef(lm3)[2]


