optimize_alpha_k_given_sigma <- function(sigma, sp, X_k, S_k, fit_km1, storage, fit_k_other_sp) {
    if(length(storage) > 0) {
        sigma_stored <- sapply(storage, "[[", "sigma")
        which_sigma_closest <- which.min(abs(log(sigma) - log(sigma_stored)))
        sigma_closest <- sigma_stored[which_sigma_closest]
        if(abs(log(sigma) - log(sigma_closest)) < 1e-6)
            return(storage[[which_sigma_closest]])
        alpha_k_init <- storage[[which_sigma_closest]]$estimate
    } else {
        if(!is.null(fit_k_other_sp)) {
            alpha_k_init <- fit_k_other_sp$alpha_k
            if(mean(abs(alpha_k_init) < 0.01))
                alpha_k_init <- rep(0.1, ncol(X_k))
        }
        else
            alpha_k_init <- rep(0.1, ncol(X_k))
    }

    fit_km1_sigma <- update_fit_sigma(fit_km1, sigma)
    opt_out <- nlm(find_pen_deviance_k, alpha_k_init,
                   sp = sp, X_k = X_k, S_k = S_k, fit_km1 = fit_km1_sigma,
                   hessian = TRUE, check.analyticals = FALSE)

    opt_out$sigma <- sigma
    opt_out
}


optimize_sigma_k <- function(sp, k, X_k, S_k, fit_km1, data, fit_k_other_sp) {
    storage <- list()

    counter <- 1
    dpen_prof_sigma <- function(sigma) {
        result <- tryCatch({
            opt_out <- optimize_alpha_k_given_sigma(sigma, sp, X_k, S_k, fit_km1, storage,
                                                    fit_k_other_sp)
            storage[[counter]] <<- opt_out
            counter <<- counter + 1
            result <- opt_out$minimum
        },  error = function(c) {
            warning(c)
            Inf
        }
        )
        result
    }


    if(!is.null(fit_k_other_sp))
        sigma_start <- fit_k_other_sp$sigma
    else
        sigma_start <- fit_km1$sigma
    
    opt_out_lsigma <- optim(sigma_start, dpen_prof_sigma, method = "L-BFGS-B",
                            lower = 1e-6, upper = fit_km1$sigma)

    opt_out <- storage[[which.min(sapply(storage, "[[", "minimum"))]]
    alpha_k <- opt_out$estimate
    sigma_hat <- opt_out$sigma

    rm(counter, storage)
        
    fit <- update_fit_sigma(fit_km1, sigma_hat)
    
    f_k <- X_k %*% alpha_k
    fit$cluster_info <- lapply(fit$cluster_info, update_cluster_info, f_k = f_k)
    fit$alpha_k <- alpha_k
    fit$l_hat <- find_loglikelihood_k(alpha_k, X_k, fit_km1, derivs = FALSE)
    fit$spr <- sp / (2  * sigma_hat^2)
    fit$lprior_fun <-  find_lprior_fun(k, alpha_k, S_k)
    fit$l_pen_hat <- - opt_out$minimum / 2
    fit$Sigma_inv <- opt_out$hessian / 2
    fit$log_ml_contrib <- approx_log_ml_contrib(fit$Sigma_inv)

    fit$f_x <- cbind(fit$f_x, f_k)
    fit$u <- find_u_hat(sigma_hat, data, fit$f0_x, fit$f_x)


    fit
}
