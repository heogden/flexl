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

    
    sp <- 10000
    sigma <- 0.1
    nbasis <- 10

    mod <- fit_given_sp(data, sp, sigma, 3, nbasis)

    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    f_0_fitted <- basis$X %*% mod[[4]]$beta_0
    f_1_fitted <- basis$X %*% mod[[4]]$beta[,1]
    f_2_fitted <- basis$X %*% mod[[4]]$beta[,2]
    f_3_fitted <- basis$X %*% mod[[4]]$beta[,3]
    sqrt(mean(f_3_fitted^2)) #' very small: so should stop at k = 2

    plot(data$x, data$y)
    points(data$x, f_0_fitted, col = 2)
    curve(mu, add = TRUE, lty = 2)

    #' fit reasonable (overfitting decreases with sp)

    x_grid <- seq(0, 1, length.out = 100)
    
    plot(data$x, f_1_fitted)
    lines(x_grid, sqrt(0.5) * delta(x_grid)[,1], lty = 2)

    plot(data$x, f_2_fitted)
    lines(x_grid, sqrt(0.1) * delta(x_grid)[,2], lwd = 5, lty = 2, col = 2)


    f_1_fun <- splinefun(data$x, f_1_fitted)
    f_2_fun <- splinefun(data$x, f_2_fitted)

    mean(f_1_fun(x_grid) * f_2_fun(x_grid))
    #' fairly close to zero: orthogonality is (approx) observed

    mean(delta(x_grid)[,1] * delta(x_grid)[,2])

    #' size of f_j functions are sensitive to choice of sigma

    mean(f_1_fun(x_grid)^2) #' close to 0.5 if we choose sigma = 0.1 (correct)
    mean(f_2_fun(x_grid)^2) #' close to 0.1 if sigma = 0.1
    
})
