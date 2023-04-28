log_det <- function(x) {
    drop_attributes(
        determinant(x, logarithm = TRUE)$modulus
    )
}


#' approximate the log marginal likelihood
#' by using a Laplace approximation
#' approximates - log phi(0, 0, Sigma)
#' if inv = TRUE, pass in the inverse of the variance (M = Sigma_inv)
#' if inv = FALSE, pass in the variance (M = Sigma)
approx_log_ml_contrib <- function(M, inv) {
    d <- nrow(M)
    ldSigma <- ifelse(inv, -1, 1) * log_det(M)
    d/2 * log(2 * pi) + 1/2 * ldSigma
}

#' approx log marginal likelihood for the last fit in the list
approx_log_ml <- function(fits) {
    fit <- fits[[length(fits)]]
    log_ml_contribs <- sapply(fits, "[[", "log_ml_contrib")
    fit$l_hat + sum(log_ml_contribs)
}
