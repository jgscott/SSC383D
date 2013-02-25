# Exercises 2
# Kernel smoothing and cross-validation. See test-kernel-smoother.R for
# example usage and testing of this code.
# Alex Reinhart

# Takes a kernel function and bandwidth `h`
weight <- function(xi, x, h, kernel) {
  kernel((xi - x) / h) / h
}

# Simple Gaussian kernel.
gaussKernel <- function(x) {
  dnorm(x)
}

# Uniform kernel:
# K(x) = 1/2 I(x), where
# I(x) = {1 for |x| <= 1, 0 otherwise}
uniformKernel <- Vectorize(function(x) {
  if (abs(x) <= 1) {
    0.5
  } else {
    0
  }
})

# Takes vectors of datapoints, a bandwidth `h`, and a chosen kernel.
# Returns a fit function. The fit function "contains" the given xs and ys,
# and predicts future values of y based on kernel smoothing of the original
# data.
# Because of Vectorize(), the fit function can take a vector of xs and return
# a vector of predicted ys.
fitKernelSmoother <- function(xs, ys, h, kernel) {
  force(xs)
  force(ys)
  Vectorize(function(x) {
    tmp <- weight(xs, x, h, kernel)
    sum(tmp * ys) / sum(tmp)
  })
}

### Cross-validation

# With the given training and testing data, fit using the given kernel at
# varying bandwidths. Return the sum of squared errors for the predictions
# vs. the testing data.
crossValidate <- function(training, testing, hs, kernel) {
  errs = c()
  for (h in hs) {
    yK <- fitKernelSmoother(training[,1], training[,2], h, kernel)
    mse <- sum((yK(testing[,1]) - testing[,2])^2) / nrow(testing)
    errs <- c(errs, mse)
  }
  cbind(hs, errs)
}

# Leave-out-one cross validation. Returns the mean squared error of the
# prediction vs. the left-out point, for all points in the dataset.
loocv <- function(data, hs, kernel) {
  errs = c()
  for (i in 1:nrow(data)) {
    e = crossValidate(data[-i,], t(as.matrix(data[i,])), hs, kernel)
    errs <- cbind(errs, e[,2])
  }
  cbind(hs, rowMeans(errs))
}
