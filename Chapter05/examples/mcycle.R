library(KernSmooth)
library(mosaic)

data(Utilities)

plot(Utilities$temp, Utilities$gasbill/Utilities$billingDays)
ll1 = locpoly(Utilities$temp, Utilities$gasbill/Utilities$billingDays, degree=0, bandwidth=5)
lines(ll1)

utilities = data.frame(temp=Utilities$temp, gasbill=Utilities$gasbill, billingdays = Utilities$billingDays)

write.csv(utilities, "utilities.csv", row.names=FALSE, quote=FALSE)

mcycle = read.table("mcycle.txt", header=TRUE)

plot(mcycle$times, mcycle$accel)

ll1 = locpoly(mcycle$times, mcycle$accel, degree=1, bandwidth=1)

plot(mcycle$times, mcycle$accel)
lines(ll1)

