log_det <- function(x) {
    drop_attributes(
        determinant(x, logarithm = TRUE)$modulus
    )
}


#' approximate the log marginal likelihood
#' by using a Laplace approximation
#' if inv = TRUE, pass in the inverse of the variance (Sigma)
#' if inv = FALSE, pass in the variance (Sigma)
approx_log_ml_contrib <- function(Sigma, inv) {
    d <- nrow(Sigma)
    ldSigma <- ifelse(inv, 1, -1) * log_det(Sigma)
    d/2 * log(2 * pi) + 1/2 * ldSigma
    
}

#' approx log marginal likelihood for the last fit in the list
approx_log_ml <- function(fits) {
    fit <- fits[[length(fits)]]
    log_ml_contribs <- sapply(fits, "[[", "log_ml_contrib")
    fit$l_hat + sum(log_ml_contribs)
}
