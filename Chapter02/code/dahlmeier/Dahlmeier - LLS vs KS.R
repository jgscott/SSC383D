### Clear previous data
rm(list = ls())

### Load Libraries

### Load Data
readdata = read.csv('utilities.csv', header = TRUE)
avgbill = readdata$gasbill / readdata$billingdays

data = cbind(readdata$temp, avgbill)
data = data[order(data[,1]),]

### Define Functions
# Normal Kernal
norm = function(x, xs, h){
  return((1/sqrt(2*pi)) * exp(-0.5*((x - xs)/h)**2))
}

# Given weights and observations, predict a new one
predict = function(w, y){
  ys = sum(w*y) / sum(w)
  return(ys)
}

# Guassian Kernal Smoothing Function
KS = function(h, x, y, xs){
#   wi = (1/h) * (1/sqrt(2*pi)) * exp(-0.5*((x - xs)/h)**2)
  wi = norm(x, xs, h)
  return(wi)
}

# Local Linear smoothing Function
LLS = function(h, x, y, xs){
  wi = norm(x, xs, h) * (sum(norm(x, xs, h) * (x - xs)**2) - (x - xs)*sum(norm(x, xs, h) * (x - xs)))
  return(wi)
}

# LOOCV
LOOCV = function(H, y, yhat){
  return(((y - yhat)/(1-H))**2)
}


### Run through to find best h's
x = data[,1]
y = data[,2]

h = seq(1, 10, .01)

yhatLLS = matrix(0, length(h), length(y))
ERRLLS = matrix(0, length(h), length(y))
LOOCVerrorLLS = NULL

yhatKS = matrix(0, length(h), length(y))
ERRKS = matrix(0, length(h), length(y))
LOOCVerrorKS = NULL

# Get LLS smoothed y values
for(i in 1:length(h)){
  for(j in 1:length(y)){
   wi = LLS(h[i],x,y,x[j])
   yhatLLS[i,j] = predict(wi, y)
   H = wi/sum(wi)
   ERRLLS[i,j] = LOOCV(H[j], y[j], yhatLLS[i,j])
  }
}

# Get KS smoothed y values
for(i in 1:length(h)){
  for(j in 1:length(y)){
    wi = KS(h[i],x,y,x[j])
    yhatKS[i,j] = predict(wi, y)
    Hi = wi/sum(wi)
    ERRKS[i,j] = LOOCV(Hi[j], y[j], yhatKS[i,j])
  }
}

# Sum the errors
LOOCVerrorLLS = apply(ERRLLS, 1, sum)
LOOCVerrorKS = apply(ERRKS, 1, sum)

# Find the best 'h' for each method
bestindexLLS = which(LOOCVerrorLLS == min(LOOCVerrorLLS))
besthLLS = h[bestindexLLS]
bestLOOCVLLS = LOOCVerrorLLS[bestindexLLS]

bestindexKS = which(LOOCVerrorKS == min(LOOCVerrorKS))
besthKS = h[bestindexKS]
bestLOOCVKS = LOOCVerrorKS[bestindexKS]

### Find confidence Intervals
var = (sd(y))**2
varLLS = NULL
varKS = NULL

for(j in 1:length(y)){
  wiLLS = LLS(h[bestindexLLS],x,y,x[j])
  varLLS[j] = var * sum(wiLLS**2) / sum(wiLLS)**2
  wiKS = KS(h[bestindexKS],x,y,x[j])
  varKS[j] = var * sum(wiKS**2) / sum(wiKS)**2
}

### Create Bounds
yU_LLS = yhatLLS[bestindexLLS,] + 1.96*sqrt(varLLS)
yL_LLS = yhatLLS[bestindexLLS,] - 1.96*sqrt(varLLS)

yU_KS = yhatKS[bestindexLLS,] + 1.96*sqrt(varKS)
yL_KS = yhatKS[bestindexLLS,] - 1.96*sqrt(varKS)

### Plot Everything
par(mfrow = c(2,2))

plot(h, LOOCVerrorLLS, main = paste("LOOCV LLS: h = ", besthLLS), xlab = 'h', ylab = 'LOOCV')
abline(v = besthLLS, col = 'red')
besthLLS
bestLOOCVLLS

plot(h, LOOCVerrorKS, main = paste("LOOCV KSS: h = ", besthKS), xlab = 'h', ylab = 'LOOCV')
abline(v = besthKS, col = 'blue')
besthKS
bestLOOCVKS

plot(x, y , main = "Data and Fit - LLS", xlab = 'x', ylab = 'f(x)')
# points(x, yhatLLS[bestindexLLS,], col = 'red')
lines(x, yhatLLS[bestindexLLS,], col = 'red')
lines(x, yU_LLS, col = 'red', lty = 2)
lines(x, yL_LLS, col = 'red', lty = 2)

plot(x, y , main = "Data and Fit - KS", xlab = 'x', ylab = 'f(x)')
# points(x, yhatKS[bestindexKS,], col = 'blue')
lines(x, yhatKS[bestindexKS,], col = 'blue')
lines(x, yU_KS, col = 'blue', lty = 2)
lines(x, yL_KS, col = 'blue', lty = 2)

par(mfrow = c(1,1))
plot(x, y, main = 'LLS, KS, and Data', xlab = 'x', ylab = 'f(x)')
lines(x, yhatLLS[bestindexLLS,], col = 'red')
lines(x, yhatKS[bestindexKS,], col = 'blue')






