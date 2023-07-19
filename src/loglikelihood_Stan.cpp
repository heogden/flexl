// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math/mix.hpp> // need this to be able to do hessian
#include <stan/math.hpp>  // pulls in everything from rev/ and prim/
#include <Rcpp.h>
#include <RcppEigen.h>


using namespace Rcpp;


// [[Rcpp::export]]
NumericVector find_loglikelihood_cluster_Stan(Eigen::VectorXd theta,
					      Eigen::MatrixXd X_c,
					      Eigen::VectorXd y_c) {
  // declarations
  double l_c;
  Eigen::VectorXd grad_l_c;
  Eigen::MatrixXd hess_l_c;

  auto n_B = X_c.cols(); // number of basis functions
  auto K = (theta.size() - 1) / n_B - 1; // number of eigenfunctions


  stan::math::hessian([X_c, y_c, n_B, K](auto theta) {
    auto sigma = theta.tail(1)[0];
    auto tau = 1 / (sigma * sigma);
    auto beta_0 = theta.head(n_B);
    auto f0_c = X_c * beta_0;
    auto z = y_c - f0_c;
    auto n = y_c.size();

    // TODO: initialize matrix with n rows and K columns
    // Or init with no columns, append them on?
    // Or could use some kind of list of the vectors, also fine.
    Matrix<stan::math::var, n, K> D;
    

    auto ldet_Sigma = -n * log(tau);
    auto Q = tau * z.dot(z);
  
    for(int k = 0; k < K; ++k) {
      auto beta_k = theta.segment(n_B * (k + 1), n_B);
      auto f_k = X_c * beta_k;
      auto b_k = tau * f_k;
      for(int j = 0; j < k; ++j) {
	auto d_j = D.col(j);
	auto d_j_f_k = d_j.dot(f_k); // Can the loop be removed and use matrix op here?
	b_k = b_k - d_j_f_k * d_j;
      }
      auto a_k = 1 + b_k.dot(f_k);
      auto d_k = b_k / stan::math::sqrt(a_k);
      D.col(k) = d_k;
  
      ldet_Sigma += log(a_k);
      auto d_k_z = d_k.dot(z);
      Q -= d_k_z * d_k_z;
    }

    return - (n * stan::math::log(2 * M_PI) + ldet_Sigma + Q) / 2;
  }, theta, l_c, grad_l_c, hess_l_c);
 
  // reformat returned result
  NumericVector l_c1 = wrap(l_c);
  NumericVector grad_l_c1 = wrap(grad_l_c);
  NumericMatrix hess_l_c1 = wrap(hess_l_c);

  l_c1.attr("gradient") = grad_l_c1;
  l_c1.attr("hessian") = hess_l_c1;
  
  return l_c1;
}
