# Exercises 3
# Simulate a Gaussian process with given covariance function.
# Explore its behavior over a range of parameters by producing animated GIFs.
# Requires the animation package.
# Alex Reinhart

library(animation)

# Return a squared exponential covariance function with given parameters
sqExpCov <- function(b, t1, t2) {
  function(x1, x2) {
    t1 * exp(-(x1 - x2)^2 / (2 * b^2)) + ifelse(x1 == x2, t2, 0)
  }
}

# For each x in xs, compute the random value f(x).
simGauss <- function(xs, b, t1, t2, r = NULL) {
  expCov = sqExpCov(b, t1, t2)

  # Build an upper triangular covariance matrix between all the xs.
  covm = matrix(nrow=length(xs), ncol=length(xs))
  for (i in 1:length(xs)) {
    for (j in i:length(xs)) {
      covm[i,j] = expCov(xs[i], xs[j])
    }
  }

  lt = lower.tri(covm)
  covm[lt] = t(covm)[lt]
  
  # Now, we simply generate a vector of standard random normals and multiply
  # by the covariance matrix's SVD.
  if (is.null(r)) {
    r = rnorm(length(xs))
  }

  s = svd(covm)

  s$u %*% diag(sqrt(s$d)) %*% t(s$u) %*% r
}

### Create animations in which we vary each parameter
xs = seq(0, 1, length.out=100)
rs = rnorm(length(xs))

ani.options(outdir = getwd(), interval = 0.1)

saveGIF({
  for (b in seq(0.01, 0.3, by=0.005)) {
    ys = simGauss(xs, b, 0.5, 10e-6, rs)
    plot(xs, ys, ylim=c(-3, 3), xlab="x", ylab="y",
         main=paste("b =", sprintf("%0.3f", b), "t1 = 0.5 t2 = 10e-6"))
  }}, movie.name="b.gif");

saveGIF({
  for (t1 in seq(0.01, 1, by=0.02)) {
    ys = simGauss(xs, 0.2, t1, 10e-6, rs)
    plot(xs, ys, ylim=c(-2, 2), xlab="x", ylab="y",
         main=paste("b = 0.2 t1 =", t1, "t2 = 10e-6"))
  }}, movie.name="t1.gif")

saveGIF({
  for (t2 in seq(0, 0.3, by=0.005)) {
    ys = simGauss(xs, 0.2, 0.5, t2, rs)
    plot(xs, ys, ylim=c(-4, 4), xlab="x", ylab="y",
         main=paste("b = 0.2 t1 = 0.5 t2=", sprintf("%0.3f", t2)))
  }}, movie.name="t2.gif")
