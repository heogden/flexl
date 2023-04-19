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

    
    sp <- 1000000
    sigma <- 0.1

    mod <- fit_given_sp(data, sp, sigma, 2, 10)

    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    f_0_fitted <- basis$X %*% mod[[3]]$beta_0
    f_1_fitted <- basis$X %*% mod[[3]]$beta[,1]
    f_2_fitted <- basis$X %*% mod[[3]]$beta[,2]

    plot(data$x, data$y)
    points(data$x, f_0_fitted, col = 2)
    curve(mu, add = TRUE, lty = 2)

    #' fit reasonable (overfitting decreases with sp)

    x_grid <- seq(0, 1, length.out = 100)
    
    plot(data$x, f_1_fitted)
    lines(x_grid, 0.7 * delta(x_grid)[,1], lty = 2)

    plot(data$x, f_2_fitted)
    lines(x_grid, 0.3 * delta(x_grid)[,2], lwd = 5, lty = 2, col = 2)

    #' constants not that important, but unsure why they take these values
    sqrt(0.5)/sqrt(0.1)
    #' is close to:
    0.68/0.3

    f_1_fun <- splinefun(data$x, f_1_fitted)
    f_2_fun <- splinefun(data$x, f_2_fitted)

    mean(f_1_fun(x_grid) * f_2_fun(x_grid))
    #' fairly close to zero: orthogonality is (approx) observed

    mean(delta(x_grid)[,1] * delta(x_grid)[,2])

    #' fitting is sensitive to choice of sigma
    #' here I have specified the correct value: but don't know
    #' this in practice!
    
})
