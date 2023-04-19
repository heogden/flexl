#' @param sp the smoothing parameter
#' @param sigma the standard deviation of the normal errors
#' @param kmax the maximum number of variation functions to use
#' @param nbasis the number of spline basis functions to use
fit_given_sp <- function(data, sp, sigma, kmax, nbasis) {
    #' find the basis to use
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    fits <- list()
    #' fit the mean-only model (k = 0)
    fits[[1]] <- fit_0(data, sp, sigma, basis)
    
    #' fit with k variation functions, fixing mean and first k-1 functions
    for(k in 1:kmax) {
        fits[[k+1]] <- fit_given_k(data, sp, sigma, k, fits[[k]], basis)
        #' stop early if f_j very close to 0
    }
    fits
}

#' @param k the number of variation functions to use 
fit_given_k <- function(data, sp, sigma, k, fit_km1, basis) {
    nbasis <- nrow(basis$S)
    
    if(k > 1) {
        #' find basis for orthogonal complement of beta_1, ..., beta_{k-1}
        transform <- find_orthogonal_complement_transform(fit_km1$beta)
            
        #' pre-compute "transformed" spline basis:
        X_k <- basis$X %*% transform
        S_k <- basis$S %*% transform
        
    } else {
        #' no restrictions needed on beta_1
        X_k <- basis$X
        S_k <- basis$S
    }
    
    #' optimize loglikelihood for alpha_k, keeping f_0, f_1, .., f_{k-1} (and sigma) fixed
    opt_out <- optim(rep(0, nbasis - k + 1), find_pen_loglikelihood_k,
                     X_k = X_k, S_k = S_k, fit_km1 = fit_km1,
                     method = "BFGS", control = list(fnscale = -1))

    alpha_k <- opt_out$par
    f_k <- X_k %*% alpha_k
    cluster_info <- lapply(fit_km1$cluster_info, update_cluster_info, f_k = f_k)

    if(k > 1) {
        beta_k <- transform %*% alpha_k
    } else {
        beta_k <- alpha_k
    }
    
    beta <- cbind(fit_km1$beta, beta_k)
    
    list(beta_0 = fit_km1$beta_0, beta = beta, cluster_info = cluster_info)
}

#' Fit the mean-only mode, and set up for use in future fits
fit_0 <- function(data, sp, sigma, basis) {
    X_0 <- basis$X
    mod_0 <- lm.fit(X_0, data$y)

    clusters <- unique(data$c)
    cluster_info <- lapply(clusters, init_cluster_info, data = data, sigma = sigma,
                           z = mod_0$residuals)
    
    list(k = 0,
         beta_0 = mod_0$coefficients,
         beta = matrix(nrow = ncol(X_0), ncol = 0),
         cluster_info = cluster_info)
}

init_cluster_info <- function(cluster, data, sigma, z) {
    rows <- which(data$c == cluster)
    n_c <- length(rows)
    Sigma_inv <- diag(1/sigma^2, nrow = n_c, ncol = n_c)
    ldet_Sigma <- 2 * n_c * log(sigma)
    z <- z[rows]
    Sigma_inv_z <- z / sigma^2
    tz_Sigma_inv_z <- sum(z^2) / sigma^2
    list(cluster = cluster, rows = rows, Sigma_inv = Sigma_inv,
         ldet_Sigma = ldet_Sigma, z = z, Sigma_inv_z = Sigma_inv_z,
         tz_Sigma_inv_z = tz_Sigma_inv_z)
}


update_cluster_info <- function(cluster_info_km1, f_k) {
    a <- f_k[cluster_info_km1$rows]
    find_info_k(a, cluster_info_km1)
}



