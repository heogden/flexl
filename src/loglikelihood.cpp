#include <Rcpp.h>

// May need to do a bit more work to process fx from R input, see
// https://teuder.github.io/rcpp4everyone_en/100_matrix.html

// Feed in subsets already, to avoid having to mess about changing indices

// [[Rcpp::export]]
double find_loglikelihood_cluster_Cpp(Rcpp::NumericVector& f0_c, Rcpp::NumericMatrix& fx_c, Rcpp::NumericVector& y_c, double sigma) {
  double tau = std::pow(sigma, -2);
  
  int n = f0_c.size();
  int K = fx_c.cols();

  Rcpp::NumericVector z = y_c - f0_c;

  Rcpp::NumericVector f_k(n);
  Rcpp::NumericVector b_k(n);

  Rcpp::NumericMatrix D(n, K);
  Rcpp::NumericVector d_j(n);
  Rcpp::NumericVector d_k(n);
  double d_j_f_k;
  double a_k;

  double ldet_Sigma = -n * log(tau);
  double Q = tau * Rcpp::sum(pow(z, 2));
  
  for(int k = 0; k < K; ++k) {
    f_k = fx_c( Rcpp::_ , k);
    b_k = tau * f_k;
    for(int j = 0; j < k; ++j) {
      d_j = D( Rcpp::_ , j);
      d_j_f_k = Rcpp::sum(d_j * f_k);
      b_k = b_k - d_j_f_k * d_j;
    }
    a_k = 1 + Rcpp::sum(b_k * f_k);
    d_k = b_k / std::sqrt(a_k);
    D( Rcpp::_ , k) = d_k;
  }
  ldet_Sigma += log(a_k);
  Q -= std::pow(Rcpp::sum(d_k * z), 2);

  return - (n * std::log(2 * M_PI) + ldet_Sigma + Q) / 2;

}
  
