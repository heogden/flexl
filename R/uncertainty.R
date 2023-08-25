find_u_sample_cluster <- function(cluster, sigma, data, f0_x, f_x) {
    inc <- which(data$c == cluster)
    f0_x_c <- f0_x[inc]
    f_x_c <- f_x[inc, , drop = FALSE]
    y_c <- data$y[inc]

    k <- ncol(f_x)
    I <- diag(1, nrow = k, ncol = k)
    A <- crossprod(f_x_c) / sigma^2 + I
    b <- crossprod(f_x_c, y_c - f0_x_c) / sigma^2
    
    u_hat <- as.numeric(solve(A, b))
    mvnfast::rmvn(1, u_hat, solve(A))
}



find_sample <- function(id, mod) {
    set.seed(id)
    par <- as.numeric(mvnfast::rmvn(1, mod$par, mod$var_par))

    par_split <- split_par(par, mod$basis$nbasis)
    
    f0_x <- mod$basis$X %*% par_split$beta0
    f0 <- find_spline_fun(par_split$beta0, mod$basis)

    if(mod$k > 0) {
        beta <- find_beta(par_split$alpha, mod$basis$nbasis, mod$k)
        f_x <- mod$basis$X %*% beta
        f <- find_spline_fun(beta, mod$basis)
    } else {
        f_x <- matrix(nrow = nrow(mod$data), ncol = 0)
        f <- NULL
    }
    
    comps <- lapply(unique(mod$data$c), find_u_sample_cluster,
                    sigma = mod$sigma, data = mod$data, f0_x = f0_x, f_x = f_x)
    u <- Reduce(rbind, comps)
    rownames(u) <- unique(data$c)
    
    list(par = par, u_hat = u, f0 = f0, f = f)
}


find_samples <- function(mod, n_samples) {
    lapply(1:n_samples, find_sample, mod = mod)
}
