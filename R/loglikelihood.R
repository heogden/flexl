find_loglikelihood_cluster <- function(rows, f0, fx, y, sigma) {
    f0_c <- f0[rows]
    fx_c <- fx[rows, , drop = FALSE]
    y_c <- y[rows]
    n_c <- length(rows)
    Sigma <- tcrossprod(fx_c, fx_c) + diag(sigma^2, nrow = n_c, ncol = n_c, names = FALSE)
    mvnfast::dmvn(y_c, mu = f0_c, sigma = Sigma, log = TRUE)
}


find_loglikelihood_beta <- function(beta0, beta, sigma, y, row_list, basis, k) {
    f0 <- as.numeric(basis$X %*% beta0)
    fx <- basis$X %*% beta
    l_contribs <- sapply(row_list, find_loglikelihood_cluster, f0 = f0, fx = fx,
                         y = y, sigma = sigma)
    sum(l_contribs)
}

split_par <- function(par, nbasis) {
    components <- rep("alpha", length(par))
    components[1:nbasis] <- "beta0"
    components[length(par)] <- "lsigma"
    split(par, components)
}

find_wiggliness_f_k <- function(beta_k, S, derivs = TRUE) {
    res <- emulator::quad.form(S, beta_k)
    if(derivs) {
        attr(res, "gradient") <- 2 * as.numeric(crossprod(beta_k, S))
        attr(res, "hessian") <-  2 * S
    }
    res
}

find_pen_deviance <- function(par, sp, y, row_list, basis, k) {
    cat("par = ", par, "\n")
    par_split <- split_par(par, basis$nbasis)

    beta0 <- par_split$beta0
    beta <- find_beta(par_split$alpha, basis$nbasis, k)$value
    sigma <- exp(par_split$lsigma)
    
    l <- find_loglikelihood_beta(beta0, beta, sigma, y, row_list, basis, k)

    beta_full <- cbind(beta0, beta)
    w_each <- apply(beta_full, 2, find_wiggliness_f_k, S = basis$S, derivs = FALSE)
    w <- sum(w_each)
    spr <- sp / (2 * sigma^2)
    pen <- - spr * w
    cat("l = ", l, ", pen = ", pen, ", dev = ", -2 * (l + pen), "\n")
  
    -2 * (l + pen)
}




