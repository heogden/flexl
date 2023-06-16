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
    
    t <- alpha_norm * (alpha_norm - alpha[1])
    dt <- 2 * (1 - alpha[1]) / alpha_norm * alpha
    dt[1] <- dt[1] - alpha_norm
    
    gamma <- 1 / t
    dgamma <- - dt / t^2

    u <- alpha
    u[1] <- u[1] - alpha_norm
    
    list(alpha = alpha,
         u = u, gamma = gamma,
         M = M)
}

#' also return derivatives wrt alpha, if requested
find_Hstarx <- function(H, x, derivs = TRUE) {
    a <- sum(Hstar$u[-1] * x)
    result <- c(0, x) - a * H_star$gamma * H_star$u
}


find_Hi <- function(alpha_i) {
    u <- alpha_i
    alpha_i_norm <- sqrt(sum(alpha_i^2))
    u[1] <- u[1] - alpha_i_norm
    new_Householder(u)
}


find_beta_i <- function(alpha_i, H_1_to_im1, P_1_to_im1) {
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
    H_list <- list()
    ##beta <- matrix(nrow = nbasis, ncol = k)
    beta <- list()
    for(i in 1:k) {
        alpha_i <- alpha[component == i]
        if(i == 1)
            beta[[i]] <- list(value = alpha_i, derivs = list(diag(nrow = length(alpha_i))))
        else
            beta[[i]] <- find_beta_i(alpha_i, H_list, P_list)
       
        H_list[[i]] <- find_Hi(alpha_i)
       
        Hi_mat <- as.matrix.Householder(H_list[[i]])
        if(i == 1)
            P_list[[i]] <-  Hi_mat
        else
            P_list[[i]] <- P_list[[i-1]] %*% Hi_mat
        
    }
    beta
}



