#include <RcppEigen.h>
#include <algorithm>    // std::max

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseVector;
typedef Map<MatrixXd>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
typedef Map<VectorXi>  MapVeci;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP sparsesgd_logit(SEXP XX, SEXP YY, SEXP NN, NumericVector pars, int npass) {
  // X is the design matrix stored in column-major format
  // i.e. with features for case i stores in column i
  // Y is the vector of number of successes per case
  // N is the vector of sample sizes, i.e. number of trials per case
  using Eigen::MappedSparseMatrix;
  using Eigen::SparseMatrix;
  const MappedSparseMatrix<double> X(as<MappedSparseMatrix<double> >(XX));
  int numobs = X.cols();
  int numfeatures = X.rows();
  const MapVeci Y(as<MapVeci>(YY));
  const MapVeci N(as<MapVeci>(NN));
  SparseVector<double> x(numfeatures);
  VectorXd beta(numfeatures);
  beta.setZero();
  NumericVector psi(numobs);
  int j,k;
  double eta = pars[0];
  double decay = pars[1];
  double gamma = eta;
  double psi0, w, delta;
  for(int pass=0; pass<npass; pass++)
    {
      k=npass*numobs;
      for(int i=0; i < numobs; i++)
	{
	  gamma = eta/pow(k+i+2,decay);
	  x = X.innerVector(i);
	  psi0 = x.dot(beta);
	  w = 1.0/(1.0+exp(-psi0));
	  psi[i] = psi0;
	  delta = Y[i] - w*N[i];
	  for (SparseVector<double>::InnerIterator it(x); it; ++it)
	    {
	      j = it.index();
	      beta(j) += gamma*delta*it.value();
	    }
	}
    }
  return List::create(Named("beta") = beta);
}
