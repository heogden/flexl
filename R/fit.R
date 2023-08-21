fit_flexl <- function(data, nbasis = 10, kmax = 10,
                      lsp_poss = -5:15) {
    if(any(is.na(data)))
        stop("There are missing values in the data, which flexl cannot handle")

    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    
    sp_poss <- exp(lsp_poss)

    fits_list <- list()
    fit_sp_poss <- list()
    
    for(i in seq_along(sp_poss)) {
        sp <- sp_poss[i]
        cat("sp = ", sp, "\n")
        if(i == 1)
            fits_other_sp <- NULL
        else
            fits_other_sp <- fits_list[[i-1]]
        
        fits <- fits_given_sp(sp, kmax, data, basis, fits_other_sp)
        if(is_k_larger_than_required(fits[[length(fits)]]))
            fit <- fits[[length(fits) - 1]]
        else
            fit <- fits[[length(fits)]]
               

        fits_list[[i]] <- fits
  
        cat("k = ", fit$k, "\n")
        fit <- add_hessian_and_log_ml(fit, basis, data)
        cat("log_ml = ", fit$log_ml, "\n")
        fit_sp_poss[[i]] <- fit

        if(i > 2)
            if(fit_sp_poss[[i]]$log_ml < fit_sp_poss[[i-1]]$log_ml)
                break         
    }

    log_ml_sp_poss <- sapply(fit_sp_poss, "[[", "log_ml")
    i_opt <- which.max(log_ml_sp_poss)
    fit_sp_poss[[i_opt]]
}

