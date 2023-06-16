find_orthogonal_complement_transform <- function(beta) {
    A <- qr(beta)
    k <- ncol(beta)
    qr.Q(A, complete = TRUE)[,-(1:k)]
}

find_beta <- function(alpha_split) {
    alpha <- alpha_split[-1]
    beta <- matrix(alpha[[1]], ncol = 1)
    if(length(alpha) > 1) {
        for(j in 2:length(alpha)) {
            transform <- find_orthogonal_complement_transform(beta)
            beta_j <- transform %*% alpha[[j]]
            beta <- cbind(beta, beta_j)
        }
    }
    beta
}

find_alpha_components <- function(nbasis, k) {
    n_each <- nbasis - 0:(k-1)
    component <- rep(1:k, times = n_each)
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
    
    a <- sum(u[-1] * x)

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


find_Hstar_x <- function(Hstar, x, deriv = TRUE) {
    a <- sum(Hstar$u[-1] * x)
    result <- c(0, x) - a * Hstar$gamma * Hstar$u
    if(deriv)
        attr(result, "gradient") <- - Hstar$gamma * outer(Hstar$u, c(0, x)) - a * Hstar$M

    result
}



find_beta_i <- function(alpha_i, Hstar_1_to_im1, P_1_to_im1) {
    step_j <- alpha_i

    i <- length(H_1_to_im1) + 1
    derivs <- P_1_to_im1[[i - 1]]
    
    for(j in rev(seq_along(H_1_to_im1))) {
        ## update: step_j <- H_j %*% (step_j)
        ## need to modify find_H0x so it also finds deriv of this wrt alpha_j
        derivs <- P_1_to_im1[[j]] %*% step_j$deriv
        step_j <- find_H0x(H_1_to_im1[[j]], step_j$value)
    }
    
    list(value = step_j$value,
         deriv = unlist(derivs))
}

find_beta <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    Hstar_list <- list()
    ##beta <- matrix(nrow = nbasis, ncol = k)
    beta <- list()
    for(i in 1:k) {
        alpha_i <- alpha[component == i]
        if(i == 1)
            beta[[i]] <- list(value = alpha_i, derivs = list(diag(nrow = length(alpha_i))))
        else
            beta[[i]] <- find_beta_i(alpha_i, Hstar_list, P_list)
       
        Hstar_list[[i]] <- find_Hstar(alpha_i)
       
        Hstar_i_mat <- as_matrix_Hstar(H_list[[i]])
        if(i == 1)
            P_list[[i]] <-  Hstar_i_mat
        else
            P_list[[i]] <- P_list[[i-1]] %*% Hstar_i_mat
        
    }
    beta
}



