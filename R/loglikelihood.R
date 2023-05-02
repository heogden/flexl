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

drop_attributes <- function(x) {
    attributes(x) <- NULL
    x
}


find_loglikelihood_k <- function(alpha_k, X_k, fit_km1, derivs = TRUE) {
    f_k <- as.numeric(X_k %*% alpha_k)
    nc <- length(fit_km1$cluster_info)

    l_comp <- c()
    if(derivs) {
        dl_comp <- matrix(NA, nrow = length(alpha_k), ncol = nc)
        d2l_comp <- array(NA, dim = c(length(alpha_k), length(alpha_k), nc))
    }
    
    for(c in seq_along(fit_km1$cluster_info)) {
        info_km1_c <- fit_km1$cluster_info[[c]]
        X_kc <- X_k[info_km1_c$rows, , drop = FALSE]
        f_kc <- f_k[info_km1_c$rows]
        lc <- ldmvnorm(f_kc, info_km1_c)
        l_comp[c] <- drop_attributes(lc)
        if(derivs) {
            dl_comp[, c] <- crossprod(X_kc, attr(lc, "gradient"))
            d2l_comp[ , , c] <- emulator::quad.form(attr(lc, "hessian"), X_kc)
        }
    }
    res <- sum(l_comp)
    if(derivs) {
        attr(res, "gradient") <- rowSums(dl_comp)
        attr(res, "hessian") <-  rowSums(d2l_comp, dims = 2)
    }
    res
}

find_wiggliness_f_k <- function(alpha_k, S_k, derivs = TRUE) {
    res <- emulator::quad.form(S_k, alpha_k)
    if(derivs) {
        attr(res, "gradient") <- 2 * as.numeric(crossprod(alpha_k, S_k))
        attr(res, "hessian") <-  2 * S_k
    }
    res
}


#' using deviance = -2 * l_pen so can minimize rather than maximise function
#' (so can use nlm for optimization)
find_pen_deviance_k <- function(alpha_k, sp, X_k, S_k, fit_km1) {
    l <- find_loglikelihood_k(alpha_k, X_k, fit_km1)
    w <-  find_wiggliness_f_k(alpha_k, S_k)
    spr <- sp / fit_km1$sigma^2
    pen <- -spr / 2 * w
    res <- -2 * (l + pen)
    attr(res, "gradient") <- - 2 * (attr(l, "gradient") - spr * attr(w, "gradient"))
    attr(res, "hessian") <- - 2 * (attr(l, "hessian") - spr * attr(w, "hessian"))
    res
}
