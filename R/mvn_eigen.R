#' Find log phi(z, 0, Sigma)
#' from eigenvalues of Sigma
ldmvnorm_eigen <- function(z, evalues) {
    n <- length(evalues)
    evalues_inv <- 1 / evalues
    const <- - n / 2 * log(2 * pi)
    log_det_Sigma_inv <- sum(log(evalues_inv))
    zT_Sigma_inv_z <- sum(z^2 * evalues_inv)
    const + log_det_Sigma_inv / 2 - zT_Sigma_inv_z / 2
}

