optimize_alpha_k_given_sigma <- function(sigma, sp, X_k, S_k, fit_km1, storage) {
    if(length(storage) > 0) {
        sigma_stored <- sapply(storage, "[[", "sigma")
        which_sigma_closest <- which.min(abs(log(sigma) - log(sigma_stored)))
        sigma_closest <- sigma_stored[which_sigma_closest]
        if(abs(log(sigma) - log(sigma_closest)) < 1e-6)
            return(storage[[which_sigma_closest]])
        alpha_k_init <- storage[[which_sigma_closest]]$par
    } else {
        alpha_k_init <- rep(0.1, ncol(X_k))
    }

    fit_km1_sigma <- update_fit_sigma(fit_km1, sigma)
    
    opt_out <- optim(alpha_k_init, find_pen_loglikelihood_k,
                     sp = sp, X_k = X_k, S_k = S_k, fit_km1 = fit_km1_sigma,
                     method = "BFGS", control = list(fnscale = -1))
    opt_out$sigma <- sigma
    opt_out
}


optimize_sigma_k <- function(sp, X_k, S_k, fit_km1) {
    storage <- list()

    counter <- 1
    lpen_prof_sigma <- function(lsigma) {
        result <- tryCatch({
            sigma <- exp(lsigma)
            opt_out <- optimize_alpha_k_given_sigma(sigma, sp, X_k, S_k, fit_km1, storage)
            storage[[counter]] <<- opt_out
            counter <<- counter + 1
            opt_out$value
        },  error = function(c) -Inf)
        result
    }
    
    opt_out_lsigma <- optim(log(fit_km1$sigma), lpen_prof_sigma, method = "BFGS",
                            control = list(fnscale = -1))

    opt_out <- storage[[which.max(sapply(storage, "[[", "value"))]]
    alpha_k <- opt_out$par
    sigma_hat <- opt_out$sigma

    rm(counter, storage)
        
    fit <- update_fit_sigma(fit_km1, sigma_hat)
    
    f_k <- X_k %*% alpha_k
    fit$cluster_info <- lapply(fit$cluster_info, update_cluster_info, f_k = f_k)
    fit$alpha_k <- alpha_k

    fit
}
