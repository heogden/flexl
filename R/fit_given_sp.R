#' @param data a data frame
#' @param sp the smoothing parameter
#' @param kmax the maximum number of variation functions to use
#' @param nbasis the number of spline basis functions to use
#' @param fve_threshold the threshold on fraction of variation explained, used to choose the number of variation functions
fit_given_sp <- function(data, sp, kmax, nbasis, fve_threshold = 1) {
    #' find the basis to use
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    fits <- list()
    #' fit the mean-only model (k = 0)
    fits[[1]] <- fit_given_k(data, sp, 0, NULL, basis)
    
    #' fit with k variation functions, fixing mean and first k-1 functions
    if(kmax > 0) {
        for(k in 1:kmax) {
            fits[[k+1]] <- fit_given_k(data, sp, k, fits[[k]], basis)
            if(find_FVE(fits)[k] > fve_threshold)
                break
        }
    }

    k <- length(fits) - 1
    fits[[k+1]]$hessian <- loglikelihood_pen_hess(fits[[k+1]]$par,
                                                  X = basis$X, y = data$y, c = data$c - 1,
                                                  sp = sp, S = basis$S, K = k)
    fits[[k+1]]$log_ml <- approx_log_ml(fits[[k+1]]$l_pen, fits[[k+1]]$hessian)
    fits
}

find_FVE <- function(fits) {
    sigmas <- sapply(fits, "[[", "sigma")
    resid_var <- min(sigmas^2)
    non_resid_var <- sigmas^2 - resid_var
    1 - non_resid_var / non_resid_var[1]
}


find_par0 <- function(fit_km1, k, nbasis) {
    if(k == 0) {
        return(c(rep(0.01, nbasis), 0))
    }
    alpha_k0 <- rep(0.01, nbasis - k + 1)

    c(fit_km1$beta0, fit_km1$alpha, alpha_k0, fit_km1$lsigma)
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
        beta <- find_beta(par_split$alpha, basis$nbasis, k)$value
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
    
    fit
}

#' fit model with k eigenfunctions, given model fit with k - 1 eigenfunctions, same sp
fit_given_fit_km1 <- function(data, sp, k, fit_km1, basis) {
    par0 <- find_par0(fit_km1, k, basis$nbasis)        
    fit_given_par0(data, sp, k, par0, basis)
   
}

#' fit model, given fit_other_sp with same k but different sp
fit_given_fit_other_sp <- function(data, sp, fit_other_sp, basis) {
    k <- fit_other_sp$k
    H <- fit_other_sp$hessian

    
    #' using hessian from fit_other_sp, do one step of Newton-Raphson to find par0
    g <- loglikelihood_pen_grad(fit_other_sp$par, basis$X, data$y, data$c - 1, sp, basis$S, k)
    par0 <- fit_other_sp$par - solve(H, g)

    fit_given_par(data, sp, k, par0, basis)
}

