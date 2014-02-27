#load libraries
library(MASS)
library(BayesBridge)

####function zone

####conditional posteriors for parameters

###conditional post. of BETA coefs. 
##inputs:
#y: n-vector, observed data points
#X: nxp design matrix, where n is the number of observations and p is the number of independent columns
#sigmasq: scalar, value for sigma squared
#tausq: scalar, value for tau squared
##Outputs:
#p-vector of beta estimates drawn from multivariate normal distribution
betacon <- function(y,X,sigmasq,tausq){
	sigma <- solve(diag(dim(X)[2])/tausq+(t(X)%*%X)/sigmasq)
	mu <- sigma%*%t(X)%*%(y/sigmasq)
	ret <- mvrnorm(1,mu=mu,Sigma=sigma)
	ret
	}

###conditional post. of tausq
##inputs:
#y: n-vector, observed data points
#BETA: p-vector, containing an estimate for each of the p columns of the design matrix
#c, d: scalars, hyperparameters on the inverse gamma prior to tau squared
##outputs:
#scalar value of tau squared drawn from inverse gamma distribution
tausqcon <- function(BETA,c,d){
	p1 <- (c+length(BETA))/2
	p2 <- (d+t(BETA)%*%BETA)/2
	ret <- 1/rgamma(1,p1,rate=p2)
	ret
	}
###conditional post. of sigmasq
##inputs:
#y: n-vector, observed data points
#X: nxp design matrix, where n is the number of observations and p is the number of independent columns
#BETA: p-vector, containing an estimate for each of the p columns of the design matrix
#a, b: scalars, hyperparameters on the inverse gamma prior to sigma squared
##outputs: 
#scalar value of sigma squared drawn from inverse gamma distribution
sigmasqcon <- function(y,X,BETA,a,b){
	p1 <- (a+length(y))/2
	error <- t(y-X%*%BETA)%*%(y-X%*%BETA)
	p2 <- (b+error)/2
	ret <- 1/rgamma(1,p1,rate=p2)
	ret
	}

###sampling function
##this function performs the gibbs sampling operation with our given conditional posteriors
##Inputs:
#nsamples: scalar, number of sampling loops performed by the gibbs sampler
#y: n-vector, observed data points
#X: nxp design matrix, where n is the number of observations and p is the number of independent columns
#sigmasq: scalar, an initial guess for the value of sigma squared
#tausq: scalar, an initial guess for the value of tau squared
#a,b,c,d: scalars, hyperparameters for the inverse gamma priors on tausq and sigmasq 
##Outputs:
#A 3-value named list:
#sigmasq: posterior draws of sigma squared
#tausq: posterior draws of tau squared
#BETA: named nsamples x p matrix containing the nsample posterior draws for the p columns of the design matrix
gibbs <- function(nsamples,y,X,sigmasq,tausq,a,b,c,d){
	ret <- list(sigmasq=mat.or.vec(nsamples,1),tausq=mat.or.vec(nsamples,1),BETA=mat.or.vec(nsamples,dim(X)[2])) # builds an output list for each parameter
	sigmatrack <- sigmasq #initial guess for sigmasq
	tautrack <- tausq #initial guess for tausq
	betatrack <- mat.or.vec(dim(X)[2],1) #create a zero vector of coefficients for BETA
	for(i in 1:nsamples){
		betatrack <- betacon(y,X,sigmatrack,tautrack) #perform a preliminary draw for the BETA coefficients
		sigmatrack <- sigmasqcon(y,X,betatrack,a,b) #use new BETA coef. draws together with data/hyperparameters to draw sigmasq
		tautrack <- tausqcon(betatrack,c,d) #use new BETA coef. draws together with data/hyperparameters to draw tausq
		ret$sigmasq[i] <- sigmatrack #assign present sigma value to "active" slot in the return list
		ret$tausq[i] <- tautrack #as above for tau
		ret$BETA[i,] <- betatrack #as above for the beta coef
	}
	if(!is.null(colnames(X))) colnames(ret$BETA) <- colnames(X) #this line makes the output much more readable by giving the BETA coefs their correct design matrix names
	ret
}

#####data analysis

#import data
data(diabetes)

#set initial hyperparameters
a <- 1
b <- 1
c <- 1
d <- 1

#set initial guesses
tauguess <- 1
sigmaguess <- 1

#build design matrix and observed output
Xd <- diabetes$x
yd <- diabetes$y
#center output values
ydc <- yd-mean(yd)

#perform gibbs sampling and check parameters for convergence
test <- gibbs(10000,ydc,Xd,sigmaguess,tauguess,a,b,c,d)

#build a standard linear model to compare against
test2 <- lm(ydc~Xd-1)

#Convergence is an issue with these sorts of algorithms, let's check
#a few plots of posterior draws against index
plot(test$sigmasq) #sigma squared convergence
plot(test$tausq) # tau squared convergence
plot(test$BETA[,1]) # estimate for contribution of age variable

#How about trending in our data? How many "independent draws" are we getting?
acf(test$sigmasq)
acf(test$tausq)
acf(test$BETA[,1])

#The posterior draws capture the sampling distribution of our estimates
#we can use simple summary statistics to compare against, for example,
#what standard linear regression would return
check1 <- apply(test$BETA,2,mean) #builds a vector denoting the means of each column of BETA
check2 <- apply(test$BETA,2,sd) #builds a vector denoting the standard deviations of each column of BETA
check <- t(rbind(check1,check2)) #binds our estimates together and makes them look like lm output
colnames(check) <- c('PostMean','PostSD') #names the new columns so we know what we're looking at
check
mean(test$sigmasq) #this should be close to the residual error from lm
mean(test$tausq) #this is our estimate of "between-variable" variance

#compare against the linear model
summary(test2)
check
sqrt(mean(test$sigmasq))
sqrt(mean(test$tausq))

#compare length of estimate vector and MSE
t(coef(test2))%*%coef(test2) #linear model
t(check[,1])%*%check[,1] #gibbs sampler

#MSE
t(resid(test2))%*%resid(test2) #linear model
t(ydc-Xd%*%check[,1])%*%(ydc-Xd%*%check[,1]) #gibbs sampler

hist(test$BETA[,5])

#fitted plot
yfit <- Xd%*%check[,1]
plot(yfit,ydc)
abline(a=0,b=1)

#residual plot
gibbsidual <- (ydc-Xd%*%check[,1])
plot(gibbsidual)

#note that with posterior draws from the distribution we can simply "look at it"
#to determine relevant summary statistics such as 95% confidence intervals

#bmi confidence interval example
bmi95 <- c(quantile(test$BETA[,3],0.025),quantile(test$BETA[,3],0.975))
bmi95
confint(test2,'Xdbmi',level=0.95)

ldl95 <- c(quantile(test$BETA[,6],0.025),quantile(test$BETA[,6],0.975))
ldl95
confint(test2,'Xdldl',level=0.95)