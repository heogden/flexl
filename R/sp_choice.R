

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
    spr <- fit$spr
    log_ml_contribs <- sapply(fits, "[[", "log_ml_contrib")
    lprior_contribs <- sapply(fits, function(x) {
        x$lprior_fun(spr)
    })
    
    fit$l_hat + sum(lprior_contribs) + sum(log_ml_contribs)
}

find_lprior_fun <- function(k, alpha_k, S_k) {
    nbasis <- length(alpha_k) + max(k - 1, 0)
    r <- min(nbasis - 2, nbasis - k + 1)
    
    #' check r is correct (can eventually remove)
    r_man <- Matrix::rankMatrix(S_k)
    if(r != r_man)
        stop("rank is wrong, r is ", r, ", r_man is ", r_man)

    w <- find_wiggliness_f_k(alpha_k, S_k, derivs = FALSE)
    ldet <- log_det_gen(S_k, r)
    lprior <- function(spr) {
        -spr * w + r/2 * log(spr) - r/2 * log(pi) + ldet / 2
    }
    
}
