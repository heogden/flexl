add_hessian_and_log_ml <- function(fit, basis, data) {
    fit$hessian <- loglikelihood_pen_hess(fit$par,
                                          X = basis$X, y = data$y, c = data$c - 1,
                                          sp = fit$sp, S = basis$S, K = fit$k)
    fit$log_ml <- approx_log_ml(fit, fit$hessian, basis)
    fit
}



split_par <- function(par, nbasis) {
    components <- rep("alpha", length(par))
    components[1:nbasis] <- "beta0"
    components[length(par)] <- "lsigma"
    split(par, components)
}

find_fit_info <- function(opt, k, basis, sp, data) {
    par <- opt$par
    l_pen <- opt$value
    par_split <- split_par(par, basis$nbasis)


    f0_x <- basis$X %*% par_split$beta0
    f0 <- find_spline_fun(par_split$beta0, basis)
    
    if(k > 0) {
        beta <- find_beta(par_split$alpha, basis$nbasis, k)
        f_x <- basis$X %*% beta
        u_hat <- find_u_hat(exp(par_split$lsigma), data, f0_x, f_x)
        f <- find_spline_fun(beta, basis)
    } else {
        f <- NULL
        u_hat <- NULL
    }

    list(k = k,
         sp = sp,
         par = par,
         l_pen = l_pen,
         opt = opt,
         beta0 = par_split$beta0,
         alpha = par_split$alpha,
         lsigma = par_split$lsigma,
         beta = beta,
         sigma = exp(par_split$lsigma),
         f0 = f0,
         f = f,
         u_hat = u_hat)
    
    
    
    
}

fit_given_par0 <- function(data, sp, k, par0, basis) {
     opt <- optim(par0, loglikelihood_pen, loglikelihood_pen_grad,
                 X = basis$X, y = data$y, c = data$c - 1,
                 sp = sp, S = basis$S, K = k,
                 method = "BFGS", control = list(fnscale = -1, maxit = 10000))
    if(opt$convergence != 0)
        warning("optim has not converged")
     fit <- find_fit_info(opt, k, basis, sp, data)
     fit <- add_hessian_and_log_ml(fit, basis, data)
    
     fit
}

#' Fit the mean-only model
fit_0 <- function(data, sp, basis) {
    X_0 <- basis$X
    
    Xt_y <- crossprod(X_0, data$y)
    XtX <- crossprod(X_0, X_0)

    
    beta_0 <- as.numeric(solve(XtX + sp * basis$S, Xt_y))
    y_hat_0 <- X_0 %*% beta_0
    resid <- data$y - y_hat_0
    sigma <- sd(resid)

    par0 <- c(beta_0, log(sigma))
    
    fit_given_par0(data, sp, 0, par0, basis)
}


#' using a fixed hessian H
optim_quasi_NR <- function(par0, H, sp, basis, k, data, grad_tol = 1e-3) {
    par_prev <- par0
    value_prev <- loglikelihood_pen(par_prev, basis$X, data$y, data$c - 1, sp, basis$S, k)

    grad_prev <- loglikelihood_pen_grad(par_prev, basis$X, data$y, data$c - 1, sp, basis$S, k)
    
    count <- 1

    convergence <- 1
    continue <- TRUE

    H_inv <- try(solve(H), TRUE)
    if(!is.matrix(H_inv))
        return(list(par = par0, value = value_prev, counts = c(count, count), convergence = 1))
    
    while(continue) {
        par <- par_prev - H_inv %*% grad_prev
        value <- loglikelihood_pen(par, basis$X, data$y, data$c - 1, sp, basis$S, k)

        grad <- loglikelihood_pen_grad(par, basis$X, data$y, data$c - 1, sp, basis$S, k)

        count <- count + 1

        if(any(!is.finite(grad))) {
            convergence <- 0
            continue <- FALSE
            par <- par_prev
            value <- value_prev
            break
        }
        
        if(mean(abs(grad)) < grad_tol) {
            convergence <- 0
            continue <- FALSE
            break
        }
        if(value <= value_prev) {
            continue <- FALSE
            par <- par_prev
            value <- value_prev
            break
        }
        if(mean(abs(grad)) > 0.9 * mean(abs(grad_prev))) {
            continue <- FALSE
            break
        }
        par_prev <- par
        value_prev <- value
        grad_prev <- grad
    }
    list(par = par, value = value, counts = c(count, count), convergence = convergence)
}





#' fit model, given fit_other_sp with same k but different sp
fit_given_fit_other_sp <- function(data, sp, fit_other_sp, basis) {
    k <- fit_other_sp$k
    H <- fit_other_sp$hessian

    
    #' using hessian from fit_other_sp, do steps of Newton-Raphson with hessian H
    opt <- optim_quasi_NR(fit_other_sp$par, H, sp, basis, k, data)

    if(opt$convergence != 0) {
        fit <- fit_given_par0(data, sp, k, opt$par, basis)
    } else {
        fit <- find_fit_info(opt, k, basis, sp, data)
    }
    
    add_hessian_and_log_ml(fit, basis, data)
}



find_par0_given_fit_km1 <- function(fit_km1, k, nbasis) {
    if(k == 0) {
        return(c(rep(0.01, nbasis), 0))
    }
    alpha_k0 <- rep(0.01, nbasis - k + 1)

    c(fit_km1$beta0, fit_km1$alpha, alpha_k0, fit_km1$lsigma)
}


#' fit model with k eigenfunctions, given model fit with k - 1 eigenfunctions, same sp
fit_given_fit_km1 <- function(data, sp, k, fit_km1, basis) {
    par0 <- find_par0_given_fit_km1(fit_km1, k, basis$nbasis)        
    fit_given_par0(data, sp, k, par0, basis)
   
}



find_par0_given_fit_kp1 <- function(fit_kp1, k, nbasis) {
    alpha_kp1_components <- find_alpha_components(nbasis, k+1)
    alpha_k <-  fit_kp1$alpha[alpha_kp1_components != (k+1)]
    
    c(fit_kp1$beta0, alpha_k, fit_kp1$lsigma)
}


#' fit model with k eigenfunctions, given model fit with k + 1 eigenfunctions, same sp
fit_given_fit_kp1 <- function(data, sp, k, fit_kp1, basis) {
    par0 <- find_par0_given_fit_kp1(fit_kp1, k, basis$nbasis)        
    fit_given_par0(data, sp, k, par0, basis)
   
}

