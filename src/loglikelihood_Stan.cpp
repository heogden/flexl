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
  auto K = (theta.size() - 1) / n_B - 1; // number of eigenfunction
  auto n = y_c.size();


  stan::math::hessian([X_c, y_c, n_B, K, n](auto theta) {
    auto sigma = theta.tail(1)[0];
    auto tau = 1 / (sigma * sigma);
    auto beta_0 = theta.head(n_B);
    std::cout << beta_0 << "\n";
    auto f0_c = X_c * beta_0;
    auto z = y_c - f0_c;

    auto ldet_Sigma = -n * log(tau);
    auto Q = tau * stan::math::dot_self(z);

    auto beta_1 = theta.segment(n_B, n_B);
    auto f_1 = X_c * beta_1;
    auto b_1 = tau * f_1;
    auto a_1 = 1 + b_1.dot(f_1);
    auto d_1 = b_1 / stan::math::sqrt(a_1);
    
    ldet_Sigma += log(a_1);
    auto d_1_z = d_1.dot(z);
    Q -= d_1_z * d_1_z;

    // std::vector<typeof(d_1)> D;
    // D.push_back(d_1);

    
    // for(int k = 1; k < K; ++k) {
    //   auto beta_k = theta.segment(n_B * (k + 1), n_B);
    //   auto f_k = X_c * beta_k;
    //   auto b_k = tau * f_k;
    //   for(int j = 0; j < k; ++j) {
    // 	auto d_j = D[j];
    // 	auto d_j_f_k = d_j.dot(f_k); 
    // 	b_k = b_k - d_j_f_k * d_j;
    //   }
    //   auto a_k = 1 + b_k.dot(f_k);
    //   auto d_k = b_k / stan::math::sqrt(a_k);
      
    //   D.push_back(d_k);
  
    //   ldet_Sigma += log(a_k);
    //   auto d_k_z = d_k.dot(z);
    //   Q -= d_k_z * d_k_z;
    // }

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
