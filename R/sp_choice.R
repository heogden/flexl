

#' approximate the log marginal likelihood
#' by using a Laplace approximation
#' approximates - log phi(0, 0, Sigma)
approx_log_ml_contrib <- function(Sigma_inv) {
    d <- nrow(Sigma_inv)
    ldSigma <- -log_det(Sigma_inv)
    d/2 * log(2 * pi) + 1/2 * ldSigma
}

#' approx log marginal likelihood for the last fit in the list
approx_log_ml <- function(fits) {
    fit <- fits[[length(fits)]]
    log_ml_contribs <- sapply(fits, "[[", "log_ml_contrib")
    lprior_contribs <- sapply(fits, "[[", "lprior_hat")
    fit$l_hat + sum(lprior_contribs) + sum(log_ml_contribs)
}

find_lprior <- function(k, alpha_k, S_k, spr) {
    #' r is nbasis - 2 
    r <- nbasis - 2
    #' check r is correct in a range of cases
    r_man <- Matrix::rankMatrix(S_k)

    if(r != r_man)
        stop("rank is wrong, r is ", r, ", r_man is ", r_man)
    
    pen <- -spr * find_wiggliness_f_k(alpha_k, S_k, derivs = FALSE)
    pen - r/2 * log(2 * pi) + 1/2 * log_det_gen(2 * spr * S_k, r)
}
