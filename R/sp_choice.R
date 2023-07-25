#' Laplace approximate log-marginal likelihood
approx_log_ml <- function(l_pen, hessian) {
    p <- nrow(hessian)
    l_pen + p/2 * log(2*pi) + 1/2 * log_det(-hessian)
}

