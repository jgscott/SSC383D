library(Rcpp)
library(microbenchmark)
library(inline)

cppFunction('
  NumericMatrix ExpKernelC(NumericVector x, NumericVector y, NumericVector pars) {
    int d1 = x.size();
    int d2 = y.size();
    double scale = pars[0];
    double range = pars[1];
    NumericMatrix C(d1,d2);
    double total = 0;
    for(int i = 0; i < d1; i++) {
    	for(int j=0; j<d2; j++)
    	{
    		C(i,j) = scale*exp(-0.5*(pow(x[i] - y[j], 2.0))/range);
    	}
    }
    return C;
  }
')


ExpKernelR = function(x,y, pars)
{
	d1 = length(x)
	d2 = length(y)
	scale = pars[1]
	range = pars[2]
	C = matrix(0, nrow=d1, ncol=d2)
	for(i in 1:d1)
	{
		for(j in 1:d2)
		{
			C[i,j] = scale*exp(-0.5*((x[i] - y[j])^2)/range)
		}
	}
	return(C)
}



x = runif(100)
y =  runif(100)
theta = c(1,1)

microbenchmark(
  ExpKernelR(x,y, theta),
  ExpKernelC(x,y, theta)
)

