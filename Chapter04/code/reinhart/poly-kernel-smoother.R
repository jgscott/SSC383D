# Exercises 2
# Local polynomial regression.
# Note: My estimates of standard deviation are wrong; I did something suspect
# in the math.
# Alex Reinhart

source("../../../Chapter02/code/reinhart/kernel-smoother.R")

sj <- function(j, xs, h, kernel) {
  function(x) {
    sum(kernel((xs - x) / h) * (xs - x)^j)
  }
}

polyWeight <- function(xi, x, s1, s2, h, kernel) {
  kernel((x - xi) / h) * (s2(x) - (xi - x) * s1(x))
}

# Fit a linear smoother to the data. Returns a fit function whose return value
# is a list with entries $val and $stds.
fitPolySmoother <- function(xs, ys, h, kernel) {
  force(xs)
  force(ys)
  s1 <- sj(1, xs, h, kernel)
  s2 <- sj(2, xs, h, kernel)
  val <- Vectorize(function(x) {
    tmp <- polyWeight(xs, x, s1, s2, h, kernel)
    sum(tmp * ys) / sum(tmp)
  })
  
  stds <- Vectorize(function(x) {
    variance <- var(xs - val(xs))
    tmp <- polyWeight(xs, x, s1, s2, h, kernel)
    sqrt(variance * sum(tmp^2) / sum(tmp)^2)
  })
  
  function(x) {
    list(val=val(x), std=stds(x))
  }
}

### Cross-validation

# With the given training and testing data, fit using the given kernel at
# varying bandwidths. Return the sum of squared errors for the predictions
# vs. the testing data.
polyCrossValidate <- function(training, testing, hs, kernel) {
  errs = c()
  for (h in hs) {
    yK <- fitPolySmoother(training[,1], training[,2], h, kernel)
    mse <- sum((yK(testing[,1]) - testing[,2])^2) / nrow(testing)
    errs <- c(errs, mse)
  }
  cbind(hs, errs)
}

# Leave-out-one cross validation. Returns the mean squared error of the
# prediction vs. the left-out point, for all points in the dataset.
polyloocv <- function(data, hs, kernel) {
  errs = c()
  for (i in 1:nrow(data)) {
    e = polyCrossValidate(data[-i,], t(as.matrix(data[i,])), hs, kernel)
    errs <- cbind(errs, e[,2])
  }
  cbind(hs, rowMeans(errs))
}
