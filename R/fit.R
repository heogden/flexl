fit_flexl <- function(data, nbasis = 10) {
    lsp_poss <- seq(-5, 10, length.out = 10)
    sp_poss <- exp(lsp_poss)

    fits_poss <- list()
    fits_poss[[1]] <- fit_given_sp(data, sp_poss[1], 10, nbasis, 0.99)
    kmax <- length(fits_poss[[1]]) - 1

    log_ml_poss <- c()
    log_ml_poss[1] <- (fits_poss[[1]])[[kmax + 1]]$log_ml
    
    for(i in 2:length(sp_poss)) {
        fits_poss[[i]] <- fit_given_sp(data, sp_poss[i], kmax, nbasis, 1, fits_poss[[i-1]])
        log_ml_poss[i] <- (fits_poss[[i]])[[kmax + 1]]$log_ml
        if(log_ml_poss[i] < log_ml_poss[i-1])
            break
    }

    #' do we need to add some larger values?
    if(i == length(sp_poss)) {
        while(log_ml_poss[i] > (log_ml_poss[i-1] + 1)) {
            lsp_poss <- c(lsp_poss, max(lsp_poss) + 1)
            sp_poss <- exp(lsp_poss)
            i <- i + 1
            fits_poss[[i]] <- fit_given_sp(data, sp_poss[i], kmax, nbasis, 1, fits_poss[[i-1]])
            log_ml_poss[i] <- (fits_poss[[i]])[[kmax + 1]]$log_ml
        }
        
    }

    log_ml_fun <- splinefun(lsp_poss, log_ml_poss)
    opt_out <- optimize(log_ml_fun, range(lsp_poss), maximum = TRUE)
    lsp <- opt_out$maximum
    sp <- exp(lsp)

    which_closest_lsp <- which.min(abs(lsp_poss - lsp))
    fit_closest_sp <- fits_poss[[which_closest_lsp]]
   
    fit <- fit_given_sp(data, sp, kmax, nbasis, 0.99, fit_closest_sp)
    FVE <- find_FVE(fit)
    index <- which(find_FVE(fit) > 0.99)[1]
    fit[[index]]                    
}

