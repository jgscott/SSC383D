library(Rcpp)
library(inline)

# A compiled C++ function, since loops in R are slow
# Inputs
# kappa: (variance, range, nugget)
#		a set of covariance parameters for the squared exponential covariance function
# x1: (n1 x 2) matrix of 2d Euclidean coordinates for spatial locations x1
# x2: (n2 x 2) matrix of 2d Euclidean coordinates for spatial locations x2
#
# Returns
# C, a matrix of dimension (n1xn2)
# C(i,j) is the covariance between f(x1[i,]) and f(x2[j,])
# where f is a Gaussian process parametrized by kappa
cppFunction('
  NumericMatrix Exp2Sigma2D(NumericMatrix x, NumericMatrix y, NumericVector kappa) {
  	double distance;
    int n1 = x.nrow();
    int n2 = y.nrow();
    int dim = x.ncol();
    NumericMatrix C(n1,n2);
    double total = 0;
    for(int i = 0; i < n1; i++) {
    	for(int j=0; j < n2; j++)
    	{
    		distance = 0;
    		for(int k=0; k<dim; k++) distance += pow(x(i,k) - y(j,k), 2.0);
    		C(i,j) = kappa[0]*exp(-0.5*distance/kappa[1]);
    		if(distance == 0.0) C(i,j) += kappa[2];
    	}
    }
    return C;
  }
')


# R versions of covariance kernels for 1d design points
# Not actually called here
expKernel <- function(x,y,kappa)
{
	return( kappa[1] * exp(-0.5*((x-y)/kappa[2])^2) + kappa[3]*(x==y) )
}

maternKernel <- function(x,y,kappa)
{
	d = sqrt(sum((x-y)^2))
	k1 = sqrt(5)*d/kappa[2]
	return( kappa[1] * exp(-k1) * (1 + k1 +5*d^2/(3*kappa[2]^2)) + kappa[3]*(x==y) )
}

expSigma <- function(kappa,t1,t2)
{
	C = outer(t1,t2,FUN = expKernel,kappa=kappa);
	return(C);
}

maternSigma <- function(kappa,t1,t2)
{
	C = outer(t1,t2,FUN = maternKernel,kappa=kappa);
	return(C);
}
