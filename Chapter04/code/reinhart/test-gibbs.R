# Exercises 4
# Test our Gibbs sampler with the diabetes dataset.
# Alex Reinhart

source("gibbs.R")
library(BayesBridge)

data(diabetes, package="BayesBridge")

X <- diabetes$x
y <- diabetes$y

# Add intercept column to design matrix.
X <- cbind(rep.int(1, nrow(X)), X)

# Compare a Gibbs fit against an ordinary linear regression model.
# Notice that most coefficients match well, but a batch of coefficients
# (tc - ltg) are way off -- possibly because they're very collinear.
gibbsFit(X, y, 10000, 1, 1, 1, 1)

lmCmp <- lm.fit(X,y)
coef(lmCmp)
