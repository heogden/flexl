find_alpha_components <- function(nbasis, k) {
    n_each <- nbasis - 0:(k-1)
    rep(1:k, times = n_each)
}


split_alpha <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    split(alpha, component)
}

find_Hstar_mat <- function(alpha) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[1] <- u[1] - alpha_norm

    t <- sum(u^2)
    
    gamma <- 2 / t

    H <- diag(nrow = length(u)) - gamma * outer(u, u)
    H[ , -1, drop = FALSE]
}


find_T_list <- function(alpha, nbasis, k) {
    alpha_list <- split_alpha(alpha, nbasis, k)

    T_list <- list()
    #' no transformation for k = 1: identity matrix
    T_list[[1]] <- diag(nrow = length(alpha_list[[1]]))
    
    if(k > 1) {
        for(j in 2:k) {
            T_list[[j]] <- T_list[[j-1]] %*% find_Hstar_mat(alpha_list[[j-1]])
        }
    }
    
    T_list
}
