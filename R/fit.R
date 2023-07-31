fit_flexl <- function(data, nbasis = 10, tFVE = 0.99) {
    if(any(is.na(data)))
        stop("There are missing values in the data, which flexl cannot handle")

    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    
    lsp_poss <- seq(-5, 10, length.out = 10)
    sp_poss <- exp(lsp_poss)

    fits_poss <- list()
    fits_poss[[1]] <- fit_given_sp_init(data, sp_poss[1], 10, basis, tFVE)

    log_ml_poss <- rep(NA, length(sp_poss))
    log_ml_poss[1] <- fits_poss[[1]]$log_ml

    
    for(i in 2:length(sp_poss)) {
        fits_poss[[i]] <- fit_given_fit_other_sp(data, sp_poss[i], fits_poss[[i-1]], basis)
        log_ml_poss[i] <- fits_poss[[i]]$log_ml       
        if(i > 2 & (log_ml_poss[i] < log_ml_poss[i-1])) {
            break
        }
        if(log_ml_poss[i] > log_ml_poss[i-1]) {
            #' consider whether we should reduce k
            #' TODO: could write a better method for refitting with all smaller k
            mod_sp <- fit_given_sp_init(data, sp_poss[i], 10, basis, fve_threshold = tFVE)
            if(mod_sp$k < fits_poss[[i]]$k) {
                fits_poss[[i]] <- mod_sp
                log_ml_poss[i] <- mod_sp$log_ml
                log_ml_poss[1:(i-1)] <- NA
            }
        }

    }

    continue <- TRUE

    #' do we need to add some larger values?
    if(i == length(sp_poss)) {
        while(continue & i < 15) {
            lsp_poss <- c(lsp_poss, max(lsp_poss) + 1)
            sp_poss <- exp(lsp_poss)
            i <- i + 1
            fits_poss[[i]] <- fit_given_fit_other_sp(data, sp_poss[i], fits_poss[[i-1]], basis)
            log_ml_poss[i] <- fits_poss[[i]]$log_ml
            #' TODO: avoid copy-paste from above: put this in separate function
            if(log_ml_poss[i] > log_ml_poss[i-1]) {
                continue <- FALSE
                mod_sp <- fit_given_sp_init(data, sp_poss[i], 10, basis, fve_threshold = tFVE)
                if(mod_sp$k < fits_poss[[i]]$k) {
                    fits_poss[[i]] <- mod_sp
                    log_ml_poss[i] <- mod_sp$log_ml
                    log_ml_poss[1:(i-1)] <- NA
                }
            }
        }
        
    }

    lsp_poss <- lsp_poss[!is.na(log_ml_poss)]
    log_ml_poss <- log_ml_poss[!is.na(log_ml_poss)]

    log_ml_fun <- splinefun(lsp_poss, log_ml_poss)
    opt_out <- optimize(log_ml_fun, range(lsp_poss), maximum = TRUE)
    lsp <- opt_out$maximum
    sp <- exp(lsp)

    id_closest <- which.min(abs(lsp_poss - lsp))

    ## TODO: should also consider whether we can reduce K
    fit_given_fit_other_sp(data, sp, fits_poss[[i]], basis)                   
}

