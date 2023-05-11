find_loglikelihood_cluster <- function(rows, f0, fx, y, sigma) {
    f0_c <- f0[rows]
    fx_c <- fx[rows, , drop = FALSE]
    y_c <- y[rows]
    n_c <- length(rows)
    Sigma <- tcrossprod(fx_c, fx_c) + diag(sigma^2, nrow = n_c, ncol = n_c, names = FALSE)
    mvnfast::dmvn(y_c, mu = f0_c, sigma = Sigma, log = TRUE)
}


find_loglikelihood <- function(alpha, sigma, y, X, row_list, nbasis, k) {
    alpha_split <- split_alpha(alpha, nbasis, k)
    beta0 <- alpha_split[[1]]
    f0 <- as.numeric(X %*% beta0)
    beta <- find_beta(alpha_split)
    fx <- X %*% beta
    l_contribs <- sapply(row_list, find_loglikelihood_cluster, f0 = f0, fx = fx,
                         y = y, sigma = sigma)
    sum(l_contribs)
}

split_alpha <- function(alpha, nbasis, k) {
    if(k == 0)
        n_each <- nbasis
    else
        n_each <- c(nbasis, nbasis - 0:(k-1))
    component <- rep(0:k, times = n_each)
    split(alpha, component)
}

find_beta <- function(alpha_split) {
    alpha <- alpha_split[-1]
    beta <- matrix(alpha[[1]], ncol = 1)
    if(length(alpha) > 1) {
        for(j in 2:length(alpha)) {
            transform <- find_orthogonal_complement_transform(beta)
            beta_j <- transform %*% alpha[[j]]
            beta <- cbind(beta, beta_j)
        }
    }
    beta
}

find_pen_deviance <- function(par, sp, y, X, row_list, nbasis, k, basis) {
    cat("par = ", par, "\n")
    sigma <- par[length(par)]
    alpha <- par[-length(par)]
    l <- find_loglikelihood(alpha, sigma, y, X, row_list, nbasis, k)
    alpha_split <- split_alpha(alpha, nbasis, k)
    beta0 <- alpha_split[[1]]
    beta <- find_beta(alpha_split)
    beta_full <- cbind(beta0, beta)
    w_each <- apply(beta_full, 2, find_wiggliness_f_k, S = basis$S, derivs = FALSE)
    w <- sum(w_each)
    spr <- sp / (2 * sigma^2)
    pen <- - spr * w
    cat("l = ", l, ", pen = ", pen, ", dev = ", -2 * (l + pen), "\n")
  
    -2 * (l + pen)
}




