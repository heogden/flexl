test_that("can fit normal model with fixed k and penalty pars", {
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

    data <- data.frame(c = c, x = x, y = y)

    sp <- 1000
    sigma <- 1

    mod <- fit_given_sp(data, sp, sigma, 1, 10)

    
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    fits <- list()
    #' fit the mean-only model (k = 0)
    fits[[1]] <- fit_0(data, sp, sigma, basis)
    fits[[2]] <- fit_given_k(data, sp, sigma, 1, fits[[1]], basis)

    fits[[2]]$beta

    f_0_fitted <- basis$X %*% fits[[2]]$beta_0
    f_1_fitted <- basis$X %*% fits[[2]]$beta

    plot(data$x, data$y)
    points(data$x, f_0_fitted, col = 2)
    curve(mu, add = TRUE, lty = 2)

    #' fit reasonable (overfitting decreases with sp)
    
    plot(data$x, f_1_fitted)

    k <- 2
    fit_km1 <- fits[[2]]
    
    nbasis <- nrow(basis$S)
    transform <- find_orthogonal_complement_transform(fit_km1$beta)
            
    #' pre-compute "transformed" spline basis:
    X_k <- basis$X %*% transform
    S_k <- basis$S %*% transform
    
})
