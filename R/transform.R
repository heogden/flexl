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
    n_each <- c(nbasis, nbasis - 0:(k-1))
    component <- rep(0:k, times = n_each)
}


split_alpha <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    split(alpha, component)
}

find_Hi <- function(alpha_i) {
    u <- alpha_i
    alpha_i_norm <- sqrt(sum(alpha_i^2))
    u[1] <- u[1] - alpha_i_norm
    new_Householder(u)
}


find_beta_i <- function(alpha_i, H_1_to_im1) {
    beta_i <- alpha_i
    for(j in rev(seq_along(H_1_to_im1)))
        beta_i <- find_H0x(H_1_to_im1[[j]], beta_i)
     
    beta_i
    
}

find_beta <- function(alpha, nbasis, k) {
    component <- find_alpha_components(nbasis, k)
    H_list <- list()
    beta <- matrix(nrow = nbasis, ncol = k)
    for(i in 1:k) {
        alpha_i <- alpha[component == i]
        if(i == 1)
            beta[,i] <- alpha_i
        else
            beta[,i] <- find_beta_i(alpha_i, H_list)
       
        H_list[[i]] <- find_Hi(alpha_i)
    }
    beta
}



