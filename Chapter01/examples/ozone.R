# Load the library
library(mlbench)

# Load the data
ozone = data(Ozone, package='mlbench')

# Look at the help file for details
?Ozone

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
y = y - mean(y)
x = as.matrix(ozone[,2:10])

# compute the estimator
betahat = solve(t(x) %*% x) %*% t(x) %*% y
