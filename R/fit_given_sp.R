#' @param sp the smoothing parameter
#' @param sigma the standard deviation of the normal errors
#' @param kmax the maximum number of variation functions to use
#' @param nbasis the number of spline basis functions to use
fit_given_sp <- function(data, sp, kmax, nbasis) {
    #' find the basis to use
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    fits <- list()
    #' fit the mean-only model (k = 0)
    fits[[1]] <- fit_0(data, sp, basis)
    fits[[1]]$log_ml <- approx_log_ml(fits)
    
    #' fit with k variation functions, fixing mean and first k-1 functions
    if(kmax > 0) {
        for(k in 1:kmax) {
            fits[[k+1]] <- fit_given_k(data, sp, k, fits[[k]], basis)
            fits[[k+1]]$log_ml <- approx_log_ml(fits)
            #' stop early if f_j very close to 0
        }
    }
    
    fits
}

#' @param k the number of variation functions to use 
fit_given_k <- function(data, sp, k, fit_km1, basis) {
    nbasis <- nrow(basis$S)
    
    if(k > 1) {
        #' find basis for orthogonal complement of beta_1, ..., beta_{k-1}
        transform <- find_orthogonal_complement_transform(fit_km1$beta)
            
        #' pre-compute "transformed" spline basis:
        X_k <- basis$X %*% transform
        S_k <- emulator::quad.form(basis$S, transform)
        
    } else {
        #' no restrictions needed on beta_1
        X_k <- basis$X
        S_k <- basis$S
    }
    
    #' optimize loglikelihood for sigma and alpha_k, keeping f_0, f_1, .., f_{k-1} (and sigma) fixed
    fit <- optimize_sigma_k(sp, X_k, S_k, fit_km1, data)

    alpha_k <- fit$alpha_k
    
    if(k > 1) {
        beta_k <- transform %*% alpha_k
    } else {
        beta_k <- alpha_k
    }
    
    fit$beta <- cbind(fit_km1$beta, beta_k)
    fit$f = find_spline_fun(fit$beta, basis)
    
    fit
}


#' Fit the mean-only mode, and set up for use in future fits
fit_0 <- function(data, sp, basis) {
    X_0 <- basis$X
    
    Xt_y <- crossprod(X_0, data$y)
    XtX <- crossprod(X_0, X_0)


    beta_0 <- as.numeric(solve(XtX + sp * basis$S, Xt_y))
    y_hat_0 <- X_0 %*% beta_0
    resid <- data$y - y_hat_0
    sigma <- sd(resid)

    clusters <- unique(data$c)
    cluster_info <- lapply(clusters, init_cluster_info, data = data, sigma = sigma,
                           z = resid)
    
    Sigma_inv <- (XtX + sp * basis$S) / sigma^2

    l_hat <- sum(dnorm(data$y, y_hat_0, sd = sigma, log = TRUE))
    #' r is rank of S
    S <- basis$S
    r <- nbasis - 2
    #' lprior_hat is the log of the prior for beta0, evaluated at beta0_hat
    spr <- sp / sigma^2
    qhat <- emulator::quad.form(S, beta_0)
    lprior_hat <- - r/2 * log(2*pi) + 1/2 * log_det_gen(spr * S, r) - spr/2 * qhat
    lpen_hat <- l_hat + lprior_hat
    
    list(k = 0,
         beta_0 = beta_0,
         beta = matrix(nrow = ncol(X_0), ncol = 0),
         sigma = sigma,
         u = matrix(nrow = length(unique(clusters)), ncol = 0),
         lpen_hat = lpen_hat,
         log_ml_contrib = approx_log_ml_contrib(Sigma_inv),
         f0_x = y_hat_0,
         f0 = find_spline_fun(beta_0, basis),
         f_x = matrix(nrow = length(y_hat_0), ncol = 0),
         f = function(x){ matrix(nrow = length(x), ncol = 0) },
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

update_fit_sigma <- function(fit_prev, sigma) {
    s2_diff <- sigma^2 - fit_prev$sigma^2

    fit_prev$sigma <- sigma
    fit_prev$cluster_info <- lapply(fit_prev$cluster_info, update_cluster_info_sigma,
                                    s2_diff = s2_diff)
    fit_prev
}


update_cluster_info_sigma <- function(cluster_info_prev, s2_diff) {
    Sigma_prev <- solve(cluster_info_prev$Sigma_inv)
    Sigma <- Sigma_prev + diag(s2_diff, nrow = nrow(Sigma_prev), ncol = ncol(Sigma_prev))
    
    Sigma_inv <- solve(Sigma)
    ldet_Sigma <- determinant(Sigma)$modulus
    
    z <- cluster_info_prev$z
    Sigma_inv_z <- as.numeric(Sigma_inv %*% z)
    tz_Sigma_inv_z <- emulator::quad.form(Sigma_inv, z)
    
    list(cluster = cluster_info_prev$cluster, rows = cluster_info_prev$rows,
         Sigma_inv = Sigma_inv, ldet_Sigma = ldet_Sigma, z = z, Sigma_inv_z = Sigma_inv_z,
         tz_Sigma_inv_z = tz_Sigma_inv_z)
    
}


