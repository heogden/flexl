test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    sp <- 100
    basis <- find_orthogonal_spline_basis(nbasis, data$x)


    k <- 3

    row_list <- split(1:nrow(data), data$c)
    alpha_components <- find_alpha_components(nbasis, k)
    
    set.seed(1)
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), 1)

    par_split <- split_par(par, basis$nbasis)

    beta0 <- par_split$beta0
    beta <- find_beta(par_split$alpha, basis$nbasis, k)$value
    sigma <- exp(par_split$lsigma)

    f0 <- as.numeric(basis$X %*% beta0)
    fx <- basis$X %*% beta
    rows <- row_list[[1]]
    y <- data$y
    
    l1 <- find_loglikelihood_cluster(rows, f0, fx, y, sigma)
    l2 <- find_loglikelihood_cluster_struc(rows, f0, fx, y, sigma)
    expect_equal(l2, l1)
    
    pen_deviance_grad_man <- numDeriv::grad(find_pen_deviance, par, sp = sp, y = data$y,
                                            row_list = row_list, basis = basis, k = k)
    pen_deviance_grad_auto <- attr(pen_deviance, "gradient")
    expect_equal(pen_deviance_grad_auto, pen_deviance_grad_man)

    pen_deviance_hess_man <- numDeriv::hessian(find_pen_deviance, par, sp = sp, y = data$y,
                                               row_list = row_list, basis = basis, k = k)
    pen_deviance_hess_auto <- attr(pen_deviance, "hessian")
    expect_equal(pen_deviance_hess_auto, pen_deviance_hess_man)
})
