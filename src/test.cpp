// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math/mix.hpp> // need this to be able to do hessian
#include <stan/math.hpp>  // pulls in everything from rev/ and prim/
#include <Rcpp.h>
#include <RcppEigen.h>


using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]


// [[Rcpp::export]]
NumericVector ldnorm_with_derivs(Eigen::VectorXd x, Eigen::VectorXd theta) {
  //declarations
  double fx;
  Eigen::VectorXd grad_fx;
  Eigen::MatrixXd hess_fx;
  
  stan::math::hessian([x](auto theta) {
    auto mu = theta[0];
    auto sigma = theta[1];
    auto contribs = stan::math::normal_lpdf(x, mu, sigma);
    return stan::math::sum(contribs);
  }, theta, fx, grad_fx, hess_fx);

  
  // reformat returned result
  NumericVector fx1 = wrap(fx);
  NumericVector grad_fx1 = wrap(grad_fx);
  NumericMatrix hess_fx1 = wrap(hess_fx);

  fx1.attr("gradient") = grad_fx1;
  fx1.attr("hessian") = hess_fx1;
  
  return fx1;
    
}
