# Exercises 2
# Examples of kernel smoothing on different kinds of data.
# Uses the smoothing functions defined in kernel-smoother.R.
# Alex Reinhart

source("kernel-smoother.R")

### Kernel smoothing
# Simulate a noisy function to test. I will use 7cos(x) + x. Then try
# different bandwidths and see how the smoothed functions compare.
xs = seq(from=0, to=30, by=0.3)
xs = sort(xs + rnorm(length(xs), sd=0.1))
ys = 7*cos(xs) + xs + rnorm(length(xs))
plot(xs, ys)

yK = fitKernelSmoother(xs, ys, 2000, gaussKernel)
lines(xs, yK(xs))

yK = fitKernelSmoother(xs, ys, 2, gaussKernel)
lines(xs, yK(xs))

yK = fitKernelSmoother(xs, ys, 1, gaussKernel)
lines(xs, yK(xs), lty=2)

yK = fitKernelSmoother(xs, ys, 0.5, gaussKernel)
lines(xs, yK(xs), lty=3, lwd=2)

### Cross-validation
# Simulate a noisy function for training. I will use 7cos(x) + x.
xs = seq(from=0, to=30, by=0.3)
xs = sort(xs + rnorm(length(xs), sd=0.2))
ys = 7*cos(xs) + xs + rnorm(length(xs))
training = cbind(xs, ys)

# For testing I'll prepare another set identically.
xs = seq(from=0, to=30, by=0.3)
xs = sort(xs + rnorm(length(xs), sd=0.2))
ys = 7*cos(xs) + xs + rnorm(length(xs))
testing = cbind(xs, ys)

crossValidate(training, testing, c(0.1, 0.5, 1, 2, 5), gaussKernel)

### Predictive validation on different datasets
xs = seq(from=0, to=1, length.out=500)

# Noisy observations of wiggly function
ys = cos(40 * xs) + 5 * xs + rnorm(length(xs), sd=0.6)
plot(xs, ys)
d = cbind(xs, ys)

# Split the sample in two random halves.
samp = sample.int(500, 250)
testing = d[samp,]
training = d[-samp,]

plot(training)
# As an example, demonstrate what happens if we try to fit every point
# with a tiny bandwidth smoother.
smooth = fitKernelSmoother(training[,1], training[,2],
                           0.002, gaussKernel)
lines(training[,1], smooth(training[,1]), col="blue")

# Now see what bandwidths do work.
b = crossValidate(training, testing, seq(from=0.1, to=1, by=0.05),
                  gaussKernel)
plot(b)

# Good observations of a wiggly function
ys = cos(40 * xs) + 5 * xs + rnorm(length(xs), sd=0.1)
d = cbind(xs, ys)
plot(xs, ys)

samp = sample.int(500, 250)
testing = d[samp,]
training = d[-samp,]

plot(crossValidate(training, testing, seq(from=0.001, to=0.06, by=0.003),
                   gaussKernel))

# Noisy observations of a smooth function
ys = cos(6 * xs) + 2 * xs + rnorm(length(xs), sd=0.4)
plot(xs, ys)
d = cbind(xs, ys)

samp = sample.int(500, 250)
testing = d[samp,]
training = d[-samp,]

crossValidate(training, testing, seq(from=0.001, to=0.5, by=0.03),
              gaussKernel)

# Good observations of a smooth function
ys = cos(6 * xs) + 2 * xs + rnorm(length(xs), sd=0.05)
plot(xs, ys)
d = cbind(xs, ys)

samp = sample.int(500, 250)
testing = d[samp,]
training = d[-samp,]

crossValidate(training, testing, seq(from=0.001, to=0.3, by=0.02),
              gaussKernel)
