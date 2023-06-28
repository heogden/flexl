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
    fits[[1]] <- fit_0(data, sp, basis)
    
    #' fit with k variation functions, fixing mean and first k-1 functions
    if(kmax > 0) {
        for(k in 1:kmax) {
            fits[[k+1]] <- fit_given_k(data, sp, k, fits[[k]], basis)
            if(find_FVE(fits)[k] > fve_threshold)
                break
        }
    }
    
    fits
}

find_FVE <- function(fits) {
    sigmas <- sapply(fits, "[[", "sigma")
    resid_var <- min(sigmas^2)
    non_resid_var <- sigmas^2 - resid_var
    1 - non_resid_var / non_resid_var[1]
}


#' @param k the number of variation functions to use 
fit_given_k <- function(data, sp, k, fit_km1, basis) {

    if(k == 0)
        opt <- fit_0(data, sp, basis)
    else {
        #' TODO: write find_par0 to add in extra zeroes
        #' for alpha_k
        par0 <- find_par0(fit_km1$par, k, basis$nbasis)

        #' TODO: redo find_pen_deviance to take these
        #' arguments only
        #' The parameter in find_pen_deviance
        #' should contain alpha0, alpha, lsigma
        opt <- nlm(par0, find_pen_deviance, sp = sp,
                   y = data$y, basis = basis)
    }
    

    #' TODO: write find_fit_info
    #' which should (eventually) include hessian,
    #' approximate log marginal likelihood
    #' stored mean and variation functions
    #' estimates of the random effects, fitted values
    #' (will need some other arguments as well as opt)
    
    find_fit_info(opt)
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

    opt <- list(estimate = c(beta_0, sigma))
    opt
}


