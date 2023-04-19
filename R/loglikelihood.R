#' Find log phi(z; 0, Sigma_k),
#' where z is fixed and stored
#' and
#' Sigma_k = Sigma_km1 + a a^T
#' for each new value of a.
#' Here info_km1 contains values needed corresponding to Sigma_km1 and z
ldmvnorm <- function(a, info_km1) {
    b <- info_km1$Sigma_inv %*% a
    c <- 1 + sum(a * b)
    ldet_Sigma <- log(c) + info_km1$ldet_Sigma
    d <- sum(info_km1$Sigma_inv_z * a)

    n <- length(info_km1$z)

    quad_form <- info_km1$tz_Sigma_inv_z - d^2/c
    -n/2 * log(2*pi) - 1/2 * ldet_Sigma - 1/2 * quad_form
}

find_info_k <- function(a, info_km1) {
    b <- as.numeric(info_km1$Sigma_inv %*% a)
    c <- 1 + sum(a * b)
    ldet_Sigma <- log(c) + info_km1$ldet_Sigma
    Sigma_inv <- info_km1$Sigma_inv - tcrossprod(b, b) / c
    z <- info_km1$z
    Sigma_inv_z <- as.numeric(Sigma_inv %*% z)
    tz_Sigma_inv_z <- sum(z * Sigma_inv_z)
    list(cluster = info_km1$cluster, rows = info_km1$rows,
         Sigma_inf = Sigma_inv, ldet_Sigma = ldet_Sigma,
         z = z, Sigma_inv_z = Sigma_inv_z, tz_Sigma_inv_z = tz_Sigma_inv_z)
    
}


find_loglikelihood_k <- function(alpha_k, X_k, fit_km1) {
    f_k <- as.numeric(X_k %*% alpha_k)
    l_comp <- c()
    for(c in seq_along(fit_km1)) {
        info_km1_c <- fit_km1$cluster_info[[c]]
        l_comp[c] <- ldmvnorm(f_k[info_km1_c$rows], info_km1_c)
    }
    sum(l_comp)
}

find_wiggliness_f_k <- function(alpha_k, S_k) {
    emulator::quad.form(S_k, alpha_k)
}

find_pen_loglikelihood_k <- function(alpha_k, X_k, S_k, fit_km1) {
    find_loglikelihood_k(alpha_k, X_k, fit_km1) - find_wiggliness_f_k(alpha_k, S_k)
}
