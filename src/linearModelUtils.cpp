#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

class Comparator {
private:
  const Rcpp::NumericVector& ref;

  bool is_na(double x) const
  {
    return Rcpp::traits::is_na<REALSXP>(x);
  }

public:
  Comparator(const Rcpp::NumericVector& ref_)
    : ref(ref_)
  {}

  bool operator()(const int ilhs, const int irhs) const
  {
    double lhs = ref[ilhs], rhs = ref[irhs];
    if (is_na(lhs)) return false;
    if (is_na(rhs)) return true;
    return lhs < rhs;
  }
};

Rcpp::NumericVector qavg_rank(Rcpp::NumericVector x)
{
  R_xlen_t sz = x.size();
  Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
  std::sort(w.begin(), w.end(), Comparator(x));

  Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      r[w[i + k]] = i + (n + 1) / 2.;
    }
  }
  return r;
}


// [[Rcpp::export]]
double aurocCPP(NumericVector score, LogicalVector boolVect, double n1,double n2) {
  NumericVector ranksOfGr = qavg_rank(score)[!boolVect];
  double U = sum(ranksOfGr) - n1 * (n1 + 1)/2;
  return 1 - U/n1/n2;
}

// [[Rcpp::export]]
NumericVector pvalLmFit(arma::vec residuals, NumericVector coefficients,int p, arma::mat qr) {
	int n = residuals.size();
	double df = n - p;

	arma::mat qr_mat = qr(span(0,p-1), span(0,p-1));
	for(int i = 0; i < p; i++) {
		for(int j = 0; j < p; j++) {
			if(i > j) qr_mat(i, j) = 0;
		}
	}

	arma::mat  qrinv = inv(qr_mat);
	double s2 = accu(pow(residuals, 2));
	s2 = s2/df;
	arma::mat vcov = qrinv * qrinv.t() * s2;
	NumericVector se = as<NumericVector>(wrap(sqrt(vcov.diag())));
	NumericVector resPT = Rcpp::pt(abs(coefficients / se), df, true, true); //Rcpp::pt( x, df, lower = true, log = false )
	for(int i = 0; i<p ;i++){
		resPT[i] =  2 * exp(log(-resPT[i]));
		resPT[i] =  Rcpp::min(NumericVector::create(resPT[i],1));
	}
	return resPT;
}

