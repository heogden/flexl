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
    

    #' increase k until optimise log ML, for smallest sp
    continue <- TRUE
    while(continue) {
        id_prev <- id
        id_k <- id_k + 1
        k <- k_poss[id_k]
        id <- ids_k_sp[id_k, id_sp]
        cat("(fit_given_fit_km1) fit: k = ", k, "sp = ", sp)
        fits_poss[[id]] <- fit_given_fit_km1(data, sp, k, fits_poss[[id_prev]], basis)
        cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")
            
        if(fits_poss[[id]]$log_ml <= fits_poss[[id_prev]]$log_ml) {
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
        cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")
        

        #' then check if we can decrease k. Continue decreasing k until
        #' we optimize log ML
        continue <- TRUE
        while(continue) {
            id_prev <- id
            id_k <- id_k - 1
            k <- k_poss[id_k]
            id <- ids_k_sp[id_k, id_sp]
            cat("(fit_given_fit_kp1) fit: k = ", k, "sp = ", sp)
            fits_poss[[id]] <- fit_given_fit_kp1(data, sp, k, fits_poss[[id_prev]], basis)
            cat(", log_ml = ", fits_poss[[id]]$log_ml, "\n")
                            
            if(fits_poss[[id]]$log_ml < fits_poss[[id_prev]]$log_ml) {
                id_k <- id_k + 1
                id <- id_prev
                continue <- FALSE
            }
            if(k == 0) {
                continue <- FALSE
            }
        }
        #' check if we have found a better log_ML at this value of sp
        #' than we had found at previous value
        if(fits_poss[[id]]$log_ml <= fits_poss[[id_prev_outer]]$log_ml) {
            continue_sp <- FALSE
            id <- id_prev_outer
        }

        if(id_sp == length(sp_poss)) {
            continue_sp <- FALSE
        }

    }
    fits_poss[[id]]                
}

