find_info_k <- function(a, info_km1) {
    b <- as.numeric(info_km1$Sigma_inv %*% a)
    c <- 1 + sum(a * b)
    ldet_Sigma <- log(c) + info_km1$ldet_Sigma
    Sigma_inv <- info_km1$Sigma_inv - tcrossprod(b, b) / c
    z <- info_km1$z
    Sigma_inv_z <- as.numeric(Sigma_inv %*% z)
    tz_Sigma_inv_z <- sum(z * Sigma_inv_z)
    list(cluster = info_km1$cluster, rows = info_km1$rows,
         Sigma_inv = Sigma_inv, ldet_Sigma = ldet_Sigma,
         z = z, Sigma_inv_z = Sigma_inv_z, tz_Sigma_inv_z = tz_Sigma_inv_z)
    
}


find_loglikelihood_k <- function(alpha_k, X_k, fit_km1) {
    f_k <- as.numeric(X_k %*% alpha_k)
    l_comp <- c()
    for(c in seq_along(fit_km1$cluster_info)) {
        info_km1_c <- fit_km1$cluster_info[[c]]
        l_comp[c] <- ldmvnorm(f_k[info_km1_c$rows], info_km1_c)
    }
    sum(l_comp)
}

find_wiggliness_f_k <- function(alpha_k, S_k) {
    emulator::quad.form(S_k, alpha_k)
}

wiggliness_f_k_grad <- function(alpha_k, S_k) {
    2 * as.numeric(crossprod(alpha_k, S_k))
}

wiggliness_f_k_hess <- function(S_k) {
    2 * S_k
}



find_pen_loglikelihood_k <- function(alpha_k, sp, X_k, S_k, fit_km1) {
    find_loglikelihood_k(alpha_k, X_k, fit_km1) - sp / fit_km1$sigma^2 * find_wiggliness_f_k(alpha_k, S_k)
}

loglikelihood_k_grad <- function(alpha_k, X_k, fit_km1) {
    f_k <- as.numeric(X_k %*% alpha_k)
    dl_comp <- matrix(NA, nrow = length(alpha_k), ncol = length(fit_km1$cluster_info))
    for(c in seq_along(fit_km1$cluster_info)) {
        info_km1_c <- fit_km1$cluster_info[[c]]
        X_kc <- X_k[info_km1_c$rows, , drop = FALSE]
        dl_comp[, c] <- crossprod(X_kc, ldmvnorm_grad(f_k[info_km1_c$rows], info_km1_c))
    }
    rowSums(dl_comp)
}

pen_loglikelihood_k_grad <- function(alpha_k, sp, X_k, S_k, fit_km1) {
    loglikelihood_k_grad(alpha_k, X_k, fit_km1) - sp / fit_km1$sigma^2 * wiggliness_f_k_grad(alpha_k, S_k)
}



loglikelihood_k_hess <- function(alpha_k, X_k, fit_km1) {
    f_k <- as.numeric(X_k %*% alpha_k)
    nc <- length(fit_km1$cluster_info)
    d2l_comp <- array(NA, dim = c(length(alpha_k), length(alpha_k), nc))
    for(c in seq_along(fit_km1$cluster_info)) {
        info_km1_c <- fit_km1$cluster_info[[c]]
        X_kc <- X_k[info_km1_c$rows, , drop = FALSE]
        d2l_comp[ , , c] <- emulator::quad.form(ldmvnorm_hess(f_k[info_km1_c$rows], info_km1_c), X_kc)
    }
    rowSums(d2l_comp, dims = 2)
}

pen_loglikelihood_k_hess <- function(alpha_k, sp, X_k, S_k, fit_km1) {
    loglikelihood_k_hess(alpha_k, X_k, fit_km1) - sp / fit_km1$sigma^2 * wiggliness_f_k_hess(S_k)
}
