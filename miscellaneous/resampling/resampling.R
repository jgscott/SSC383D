library(mosaic)

# Nonparametric bootstrap
heights = read.csv("heights.csv")

plot(SHGT ~ MHGT, data=heights)

lmM = lm(SHGT ~ MHGT, data=heights)
lmF = lm(SHGT ~ FHGT, data=heights)

N = nrow(heights)

NMC = 500
betasave = matrix(0, NMC, 2)
for(i in 1:NMC) {
	keep = sample(1:N, N, replace=TRUE)
	lmstar = lm(SHGT ~ MHGT, data=heights[keep,])
	betasave[i,] = coef(lmstar)
}

hist(betasave[,1])
hist(betasave[,2])


########
# Residual-resampling bootstrap
########

chym = read.csv("chymotrypsin.csv")
plot(Rate ~ Conc, data=chym)

# Michaelis-Menten kinetics

mmpredict = function(x, vmax, km) {
	vmax*x/{km+x}
}

target = function(theta, x, y) {
	vmax = exp(theta[1])
	km = exp(theta[2])
	ypred = mmpredict(x, vmax, km)
	sum(0.5*{y-ypred}^2)
}

mymax = optim(c(0,0), target, x = chym$Conc, y = chym$Rate)
theta = mymax$par

rgrid = seq(0,0.5, length=100)
plot(Rate ~ Conc, data=chym)
lines(rgrid, mmpredict(rgrid, exp(theta[1]), exp(theta[2])))

### Now bootstrap
yhat = mmpredict(chym$Conc, exp(theta[1]), exp(theta[2]))
eps = chym$Rate - yhat
N = nrow(chym)

NMC = 250
thetasave = matrix(0, NMC, 2)

for(i in 1:NMC) {
	estar = sample(eps, N, replace=TRUE)
	ystar = yhat + estar
	mymax = optim(c(0,0), target, x = chym$Conc, y = ystar)
	thetastar = mymax$par
	plot(chym$Conc, ystar,ylim=c(1,2.1))
	lines(rgrid, mmpredict(rgrid, exp(thetastar[1]), exp(thetastar[2])))
	thetasave[i,] = theta
}

# Inspect the sampling distributions
hist(exp(thetasave[,1]))
hist(exp(thetasave[,2]))
plot(Rate ~ Conc, data=chym)
lines(rgrid, mmpredict(rgrid, exp(theta[1]), exp(theta[2])))




## Parametric bootstrap
brca = read.csv('brca.csv', header=TRUE)

N = nrow(brca)

glm1 = glm(recall ~ factor(radiologist) + age5059 + age6069 + age70plus + familyhistory + biopsurg + symptoms + premeno + postmenohormone + postmenounknown + previousmammogram + density1 + density3 + density4, data=brca, family=binomial)


NMC = 500
betasave = matrix(0, NMC, length(coef(glm1)))
wpred = fitted(glm1)
for(i in 1:NMC) {
	ystar = rbinom(N, 1, wpred)
	glm.boot = glm(ystar ~ factor(radiologist) + age5059 + age6069 + age70plus + familyhistory + biopsurg
		+ symptoms + premeno + postmenohormone + postmenounknown + previousmammogram + density1 + density3
		+ density4, data=brca, family=binomial)
	betasave[i,] = coef(glm.boot)
}
colnames(betasave) = names(coef(glm1))
hist(betasave[,1])
hist(betasave[,2])

boot.se = apply(betasave, 2, sd)
asymp.se = sqrt(diag(vcov(glm1)))

cbind(boot.se, asymp.se)

hist(betasave[,16],50)




####

gardasil = read.csv('gardasil.csv')

glm1 = glm(Completed ~ AgeGroup + InsuranceType, data=gardasil, family='binomial')
glm2 = glm(Completed ~ AgeGroup + InsuranceType + Location, data=gardasil, family='binomial')


