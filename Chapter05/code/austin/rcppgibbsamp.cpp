/*A few notes concerning the manner in which Rcpp operates:
1. the Rcpp::depends block at the top tells Rcpp that it needs to use the RcppArmadillo package during compiling or it'll screw up
2. *every* function that you want to export to R must be prefaced with the noted Rcpp::export comment block written exactly as seen for the gibbsC function. Again, you must put this comment in front of each function that you want exported to R or it will be ignored.
3. mat and vec are matrix and column vector respectively, 
4. read the armadillo API documentation if you need to do linear algebra problems
5. the wrap() function converts, from what I can tell, every object given to it to SEXP, rendering RCPP capable of returning the output to R in most cases.
 */

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace arma; 
using namespace Rcpp;

vec mvrnormArma(vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  vec Y = randn<vec>(ncols);
  mat eigvec;
  vec eigval;
  eig_sym(eigval,eigvec,sigma);
  mat L = diagmat(sqrt(eigval));
  L = eigvec*L;
  return L*Y+mu;
}

vec betaconC(vec y, mat X, double sigmasq, double tausq){
  mat ident;
  ident.eye(X.n_cols,X.n_cols);
  mat sigma = inv(ident/tausq+(X.t()*X)/sigmasq);
  vec mu = sigma*X.t()*(y/sigmasq);
  vec ret = mvrnormArma(mu,sigma);
  return(ret);
}

double tausqconC(vec BETA, double c, double d){
  RNGScope scope;
  double p1 = as_scalar((c+BETA.n_rows)/2.0);
  double p2 = as_scalar((d+BETA.t()*BETA)/2.0);
  double gam = R::rgamma(p1,1.0/p2); // big note, gamma is parameterized by scale
  double ret = 1.0/gam;
  return(ret);
}


double sigmasqconC(vec y, mat X, vec BETA, double a, double b){
  RNGScope scope;
  double p1 = as_scalar((a+y.n_rows)/2.0);
  double error = as_scalar(trans(y-X*BETA)*(y-X*BETA));
  double p2 = (b+error)/2.0;
  double gam = R::rgamma(p1,1.0/p2);
  double ret = 1.0/gam;
  return(ret);
}

// [[Rcpp::export]]

SEXP gibbsC(int nsamples, vec y, mat X, double sigmasqguess, double tausqguess,
            double a, double b, double c, double d){
  vec sigmareturn;
  sigmareturn.zeros(nsamples);
  vec taureturn;
  taureturn.zeros(nsamples);
  mat betareturn;
  betareturn.zeros(nsamples,X.n_cols);

  double sigmatrack = sigmasqguess;
  double tautrack = tausqguess;
  vec betatrack;
  betatrack.zeros(X.n_cols);

  for(int i=0; i<nsamples ; i++){
    betatrack = betaconC(y,X,sigmatrack,tautrack);
    sigmatrack = sigmasqconC(y,X,betatrack,a,b);
    tautrack = tausqconC(betatrack,a,b);
    betareturn.row(i) = trans(betatrack);
    sigmareturn(i) = sigmatrack;
    taureturn(i) = tautrack;
  } 
  List ret;
  ret["sigmasq"] = wrap(sigmareturn);
  ret["tausq"] = wrap(taureturn);
  ret["BETA"] = wrap(betareturn); 
  return(ret);
}
