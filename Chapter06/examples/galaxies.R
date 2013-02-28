# Simulate some test scores

N = 650
w = c(100, 400, 150); w = w/sum(w)
gamma = sample(1:3, size=N, prob=w, replace=TRUE)
theta = c(55, 70, 85)
tau = c(15, 10, 5)

y = rnorm(N, theta[gamma], tau[gamma])
hist(y,30)

# Galaxy data
galaxies = scan("galaxies.txt")
