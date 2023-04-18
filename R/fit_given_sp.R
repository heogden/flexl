#' @param sp the smoothing parameter
#' @param kmax the maximum number of variation functions to use
#' @param nbasis the number of spline basis functions to use
fit_given_sp <- function(data, sp, kmax, nbasis) {
    #' find the basis to use
    
    #' fit the mean-only model (k = 0)

    #' fit with k variation functions, fixing mean and first k-1 functions
    for(k in 1:kmax) {
        #' stop early if f_j very close to 0
    }
    #' (can extract fit with k variation functions by considering subset of
    #'  whole parameter vector for last fit)
}

#' @param beta_lk k x nbasis matrix with spline coefficients for each of f_0, .., f_{k-1}
fit_given_k <- function(data, k, beta_lk, basis) {
    #' k = 0 is mean only
    if(k > 1) {
        #' find basis for orthogonal complement of beta_1, ..., beta_km1
        

    } else {
        #' no restrictions needed on beta_0 and beta_1
    }

    #' pre-compute "transformed" spline basis:
    #'   A = spline design (full data) %*% transform mat
    #'   and S_transform = t(transform mat) %*% S %*% transform mat
    
    #' optimize loglikelihood for alpha_k, keeping beta_mk fixed

    #' within loglikelihood:
    #'     f_k(t) = A %*% alpha_k (for full data set)
    #' within penalty:
    #'   only need to compute expected wiggliness w_k of f_k
    #'     (others w_0, ..., w_{k-1} already stored and are constants)
    #'   w_k = t(alpha_k) %*% S_transform %*% alpha_k

    #' to calculate loglikelihood, then need to split by cluster.
    #' For each cluster, pass in information from fit_km1
    #' This information will be a list, with one element for each cluster
    #' So need to pass in relevant part of that list, along with subsets
    #' of z and f_k for each cluster
    
    #' As well as returning optimal beta_k, should return other quantities
    #' which will be needed for next k. Group beta_k and all this extra information
    #' as fit_k, to be passed on. So actually should pass in fit_km1, rather than beta_lm

    #' Return list of information, with one element for each cluster, including
    #' (all at optimal beta_k):
    #'   Sigma_k_inv
    #'   log(det(Sigma_k))
    #'   z^T Sigma_k_inv z
    #'   Sigma_k_inv z


    #' Globally (not separately for each cluster), also need
    #'   expected wiggliness of f_0, ..., f_k
    #'   z = y - mu (can be found after k = 0, once we know mu, then passed on)


}

