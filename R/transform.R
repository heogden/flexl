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
    if(k == 0)
        n_each <- nbasis
    else
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


find_beta <- function(alpha, nbasis, k) {
    
}
