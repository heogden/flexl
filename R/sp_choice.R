correct_lprior_alpha <- function(fit, basis) {
    k <- fit$k
    if(k > 1) {
        T_list <- find_T_list(fit$alpha, basis$nbasis, fit$k)
    }
    
    S <- basis$S
    contrib <- c()

    for(j in 0:k) {
        if(j > 1) {
            T_j <- T_list[[j]]
            S_j <- t(T_j) %*% S %*% T_j
        } else {
            S_j <- S
        }
        r_j <- min(nbasis - 2, nbasis - j + 1)
        ldg_S_j_spr <- log_det_gen(S_j, r_j) + r_j * (log(fit$sp) - 2 * fit$lsigma)
        contrib[j+1] <- -r_j * log(2 * pi) / 2 + ldg_S_j_spr / 2
    }

    sum(contrib)
}

#' Laplace approximate log-marginal likelihood
approx_log_ml <- function(fit, hessian, basis) {
    p <- nrow(hessian)
    cat("terms are: \n")
    cat("l_pen = ", fit$l_pen, "\n")
    cat("correction = ", correct_lprior_alpha(fit, basis), "\n")
    cat("p/2 * log(2*pi) = ",  p/2 * log(2*pi), "\n")
    cat("1/2 * log_det(hessian) = ", 1/2 * log_det(hessian), "\n")
    fit$l_pen + correct_lprior_alpha(fit, basis) + p/2 * log(2*pi) - 1/2 * log_det(hessian)
}

