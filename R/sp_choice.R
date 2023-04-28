

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
    fit$lpen_hat + sum(log_ml_contribs)
}
