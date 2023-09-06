find_u_hat_cluster <- function(cluster, sigma, data, f0_x, f_x) {
    inc <- which(data$c == cluster)
    f0_x_c <- f0_x[inc]
    f_x_c <- f_x[inc, , drop = FALSE]
    y_c <- data$y[inc]

    k <- ncol(f_x)
    I <- diag(1, nrow = k, ncol = k)
    A <- crossprod(f_x_c) / sigma^2 + I
    b <- crossprod(f_x_c, y_c - f0_x_c) / sigma^2
    as.numeric(solve(A, b))
}


find_u_hat <- function(sigma, data, f0_x, f_x) {
    clusters <- sort(unique(data$c))
    comps <- lapply(clusters, find_u_hat_cluster,
                    sigma = sigma, data = data, f0_x = f0_x, f_x = f_x)
    u_hat <- Reduce(rbind, comps)
    rownames(u_hat) <- clusters
    u_hat
}
