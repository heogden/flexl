find_alpha_components <- function(nbasis, k) {
    n_each <- nbasis - 0:(k-1)
    rep(1:k, times = n_each)
}


split_alpha <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    split(alpha, component)
}

find_Hstar <- function(alpha) {
    alpha_norm <- sqrt(sum(alpha^2))

    u <- alpha
    u[1] <- u[1] - alpha_norm

    t <- sum(u^2)
    dt <- 4 * alpha - 2 * alpha[1] * alpha / alpha_norm
    dt[1] <- dt[1] - 2 * alpha_norm
    
    gamma <- 2 / t
    dgamma <- - 2 * dt / t^2
    
    n <- length(alpha)
    du <- diag(nrow = n)
    du[1,] <- du[1,] - alpha / alpha_norm

    M <- outer(u, dgamma) + gamma * du
    
    list(u = u, gamma = gamma,
         M = M)
}

as_matrix_Hstar <- function(Hstar) {
    u <- Hstar$u
    gamma <- Hstar$gamma

    H <- diag(nrow = length(u)) - gamma * outer(u, u)
    H[ , -1, drop = FALSE]
}


find_Hstar_x <- function(Hstar, x) {
    a <- sum(Hstar$u[-1] * x)
    list(value = c(0, x) - a * Hstar$gamma * Hstar$u,
         gradient = - Hstar$gamma * outer(Hstar$u, c(0, x)) - a * Hstar$M)
}



find_beta_i <- function(alpha_i, Hstar_list, P_list) {

    i <- length(Hstar_list) + 1

    gradients <- list()
    gradients[[i]] <- P_list[[i - 1]]

    s_j <- list(value = alpha_i)
    
    for(j in (i-1):1) {
        s_j <- find_Hstar_x(Hstar_list[[j]], s_j$value)
        if(j == 1)
            gradients[[j]] <- s_j$gradient
        else
            gradients[[j]] <- P_list[[j - 1]] %*% s_j$gradient
    }
    
    list(value = s_j$value,
         gradient = Reduce(cbind, gradients))
}

find_beta <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    Hstar_list <- list()
    P_list <- list()
    beta <- matrix(nrow = nbasis, ncol = k)
    gradient <- array(0, dim = c(nbasis, length(component), k)) 
        
    for(i in 1:k) {
        alpha_i <- alpha[component == i]
        if(i == 1)
            beta_i <- list(value = alpha_i, gradient = diag(nrow = length(alpha_i)))
        else
            beta_i <- find_beta_i(alpha_i, Hstar_list, P_list)

        beta[,i] <- beta_i$value
        gradient[, (component <= i), i] <- beta_i$gradient
        
        Hstar_list[[i]] <- find_Hstar(alpha_i)
       
        Hstar_i_mat <- as_matrix_Hstar(Hstar_list[[i]])
        if(i == 1)
            P_list[[i]] <-  Hstar_i_mat
        else
            P_list[[i]] <- P_list[[i-1]] %*% Hstar_i_mat
        
    }
    list(value = beta,
         gradient = gradient)
}



