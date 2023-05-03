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
    alpha <- as.numeric(t(sapply(lambda, function(lambda_i){rnorm(d, sd = sqrt(lambda_i))})))
    eta <-  m + as.numeric(Z %*% alpha)

    sigma <- 0.1
    y <- eta + rnorm(n, sd = sigma)

    list(data = data.frame(c = c, x = x, y = y),
         mu = mu,
         delta = delta,
         eta = eta)
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
    alpha <- as.numeric(t(sapply(lambda, function(lambda_i){rnorm(d, sd = sqrt(lambda_i))})))
    eta <-  m + as.numeric(Z %*% alpha)

    sigma <- 0.1
    y <- eta + rnorm(n, sd = sigma)

    list(data = data.frame(c = c, x = x, y = y),
         mu = mu,
         delta = delta,
         eta = eta)
}

