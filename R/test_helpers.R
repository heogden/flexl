generate_test_data_1 <- function() {
    mu <- function(x) {
        x
    }
    delta <- function(x) {
        delta_1 <- rep(1, length(x))
        delta_2 <- 2 * sqrt(3) * (x - 1/2)
        unname(cbind(delta_1, delta_2))
    }

    set.seed(1)
    d <- 50
    n_c <- 10
    n <- d * n_c
    c <- rep(1:d, each = n_c) 
    x <- runif(n)

    m <- mu(x)


    Z_list <- tapply(x, c, delta)
    Z <- Matrix::bdiag(Z_list)
    
    lambda <- c(0.5, 0.1)
    K <- length(lambda)
    alpha_mat <- t(sapply(lambda, function(lambda_i){rnorm(d, sd = sqrt(lambda_i))}))
    alpha <- as.numeric(alpha_mat)
    eta <-  m + as.numeric(Z %*% alpha)

    eta_fun <- function(x, c) { 
        Z <- delta(x)
        mu(x) + as.numeric(Z %*% alpha_mat[,c])
    }
    


    sigma <- 0.1
    y <- eta + rnorm(n, sd = sigma)

    list(data = data.frame(c = c, x = x, y = y),
         mu = mu,
         delta = delta,
         eta = eta,
         eta_fun = eta_fun)
}

generate_test_data_2 <- function() {
    mu <- function(x) {
        sin(2 * pi * x)
    }
    delta <- function(x) {
        delta_1 <- rep(1, length(x))
        delta_2 <- 2 * sqrt(3) * (x - 1/2)^2
        unname(cbind(delta_1, delta_2))
    }

    set.seed(1)
    d <- 50
    n_c <- 10
    n <- d * n_c
    c <- rep(1:d, each = n_c) 
    x <- runif(n)

    m <- mu(x)


    Z_list <- tapply(x, c, delta)
    Z <- Matrix::bdiag(Z_list)
    
    lambda <- c(0.5, 0.1)
    K <- length(lambda)
    alpha_mat <- t(sapply(lambda, function(lambda_i){rnorm(d, sd = sqrt(lambda_i))}))
    alpha <- as.numeric(alpha_mat)
    eta <-  m + as.numeric(Z %*% alpha)

    eta_fun <- function(x, c) { 
        Z <- delta(x)
    
        mu(x) + as.numeric(Z %*% alpha_mat[,c])
    }
    

    sigma <- 0.1
    y <- eta + rnorm(n, sd = sigma)

    list(data = data.frame(c = c, x = x, y = y),
         mu = mu,
         delta = delta,
         eta = eta,
         eta_fun = eta_fun)
}

generate_test_data_0 <- function() {
    mu <- function(x) {
        sin(2 * pi * x)
    }
    
    set.seed(1)
    d <- 50
    n_c <- 10
    n <- d * n_c
    c <- rep(1:d, each = n_c) 
    x <- runif(n)

    m <- mu(x)
    
    sigma <- 0.1
    y <- m + rnorm(n, sd = sigma)

    list(data = data.frame(c = c, x = x, y = y),
         mu = mu)
}


#' @export
simulate_ri <- function(seed, beta0, beta1, sigma_u, sigma, n_clusters, n_obs_per_cluster) {
    set.seed(seed)
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    x <- runif(length(c))

    u <- rnorm(n_clusters, sd = sigma_u)

    mu_c <- function(x, c) {
        beta0 + beta1 * x + u[c]
    }

    
    pred_data <- expand.grid(x = seq(min(x), max(x), length.out = 100),
                             c = 1:n_clusters)
    pred_data$mu_c <- mu_c(pred_data$x, pred_data$c)
   
    mu <- beta0 + beta1 * x + u[c]
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon

    data <- data.frame(c = c,
                       x = x,
                       y = y,
                       mu = mu)
    list(data = data, pred_data = pred_data)
}

#' @export
simulate_rs <- function(seed, beta0, beta1, sigma_u, sigma_u_slope, corr_ri_rs, sigma, n_clusters, n_obs_per_cluster) {
    set.seed(seed)
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    x <- runif(length(c))

    Sigma_u_12 <- corr_ri_rs * sigma_u * sigma_u_slope
    Sigma_u <- matrix(c(sigma_u^2, Sigma_u_12, Sigma_u_12, sigma_u_slope^2), nrow = 2, ncol = 2)
    
    u <- mvtnorm::rmvnorm(n_clusters, sigma = Sigma_u)

    mu_c <- function(x, c) {
        (beta0 + u[c, 1]) + (beta1 + u[c, 2]) * x
    }

    
    pred_data <- expand.grid(x = seq(min(x), max(x), length.out = 100),
                             c = 1:n_clusters)
    pred_data$mu_c <- mu_c(pred_data$x, pred_data$c)
   
    mu <- (beta0 + u[c, 1]) + (beta1 + u[c, 2]) * x
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon

    data <- data.frame(c = c,
                       x = x,
                       y = y,
                       mu = mu)
    
    list(data = data, pred_data = pred_data)
}

#' @export
simulate_1dv <- function(seed, beta0, beta1, sigma_u, sigma, n_clusters, n_obs_per_cluster) {
    set.seed(seed)
    c <- rep(1:n_clusters, each = n_obs_per_cluster)
    x <- runif(length(c))

    u <- rnorm(n_clusters, sd = sigma_u)

    s <- function(x) {
        1/(0.5 + 0.1 * (10 * x - 5)^2)
    }
    
    mu_c <- function(x, c) {
        s(x) * (beta0 + beta1 * x + u[c])
    }

    
    pred_data <- expand.grid(x = seq(min(x), max(x), length.out = 100),
                             c = 1:n_clusters)
    pred_data$mu_c <- mu_c(pred_data$x, pred_data$c)
   
    mu <- s(x) * (beta0 + beta1 * x + u[c])
    epsilon <- rnorm(length(mu), sd = sigma)
    
    y <- mu + epsilon

    data <- data.frame(c = c,
                       x = x,
                       y = y,
                       mu = mu)
    list(data = data, pred_data = pred_data)
}
