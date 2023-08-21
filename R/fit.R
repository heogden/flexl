fit_flexl <- function(data, nbasis = 10, kmax = 10,
                      lsp_poss = -5:15) {
    if(any(is.na(data)))
        stop("There are missing values in the data, which flexl cannot handle")

    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    
    sp_poss <- exp(lsp_poss)
    k_poss <- 0:kmax

    ids_k_sp <- matrix(1:(length(sp_poss) * length(k_poss)),
                       nrow = length(k_poss),
                       ncol = length(sp_poss))

    fits_poss <- list()
        
    #' start with first sp and k = 0
    id_sp <- 1
    sp <- sp_poss[id_sp]
    id_k <- 1
    id <- ids_k_sp[id_k, id_sp]

    cat("(fit_0) fit: k = 0, sp = ", sp)
    fits_poss[[id]] <- fit_0(data, sp, basis)
    cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")
    

    
    #' increase k until no change in sigma, for smallest sp
    continue <- TRUE
    while(continue) {
        id_prev <- id
        id_k <- id_k + 1
        k <- k_poss[id_k]
        id <- ids_k_sp[id_k, id_sp]
        cat("(fit_given_fit_km1) fit: k = ", k, "sp = ", sp)
        fits_poss[[id]] <- fit_given_fit_km1(data, sp, k, fits_poss[[id_prev]], basis)
        cat(", sigma = ", fits_poss[[id]]$sigma)
        cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")
        
        tol <- 1e-5
        if(fits_poss[[id]]$sigma/fits_poss[[id_prev]]$sigma > (1 - tol)) {
            id_k <- id_k - 1
            id <- id_prev
            continue <- FALSE
        } 
        
        if(k == kmax) {
            continue <- FALSE
        }
    }

    continue_sp <- TRUE
    while(continue_sp) {
        #' increase sp to next value
        id_prev_outer <- id
    
        id_sp <- id_sp + 1
        sp <- sp_poss[id_sp]
        
        id <- ids_k_sp[id_k, id_sp]
        k <- k_poss[id_k]
        cat("(fit_given_fit_other_sp) fit: k = ", k, "sp = ", sp)
        fits_poss[[id]] <- fit_given_fit_other_sp(data, sp, fits_poss[[id_prev_outer]], basis)
        cat(", sigma = ", fits_poss[[id]]$sigma)
        cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")
        

        #' then check if we can decrease k, without making sigma smaller
        continue <- TRUE
        while(continue) {
            id_prev <- id
            id_k <- id_k - 1
            k <- k_poss[id_k]
            id <- ids_k_sp[id_k, id_sp]
            cat("(fit_given_fit_kp1) fit: k = ", k, "sp = ", sp)
            fits_poss[[id]] <- fit_given_fit_kp1(data, sp, k, fits_poss[[id_prev]], basis)
            cat(", sigma = ", fits_poss[[id]]$sigma)
            cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")


            
            if(fits_poss[[id_prev]]$sigma/fits_poss[[id]]$sigma < (1 - tol)) {
                id_k <- id_k + 1
                id <- id_prev
                continue <- FALSE
            }
            if(k == 0) {
                continue <- FALSE
            }
        }

  
        
        #' check if we have found a better log_ML at this value of sp
        #' than we had found at previous value of sp (for same k)
        id_prev <- ids_k_sp[id_k, id_sp - 1]
        if(is.null(fits_poss[[id_prev]])) {
            sp_prev <- sp_poss[id_sp - 1]
            fits_poss[[id_prev]] <- fit_given_fit_other_sp(data, sp_prev, fits_poss[[id]], basis)
        }
        
        
        #if(fits_poss[[id]]$log_ml <= fits_poss[[id_prev]]$log_ml) {
        #    continue_sp <- FALSE
        #    id <- id_prev
        #}

        if(id_sp == length(sp_poss)) {
            continue_sp <- FALSE
        }

    }
    log_ml_poss <- sapply(fits_poss, "[[", "log_ml")
    log_ml_poss[sapply(log_ml_poss, is.null)] <- -Inf
    fits_poss[[which.max(log_ml_poss)]]
}

