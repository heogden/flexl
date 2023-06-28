test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    sp <- 100
    basis <- find_orthogonal_spline_basis(nbasis, data$x)


    k <- 2

    row_list <- split(1:nrow(data), data$c)
    alpha_components <- find_alpha_components(nbasis, k)
    
    set.seed(1)
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), 0)
    
    pen_deviance <- find_pen_deviance(par, sp, data$y, row_list, basis, k)

    opt_out <- optim(par, find_pen_deviance, sp = sp, y = data$y, row_list = row_list,
                     basis = basis, k = k, method = "L-BFGS-B", lower = -5, upper = 5,
                     control = list(maxit = 1000))

    #' check: is the fit reasonable?
    #' this is easiest to do if we can do prediction. So next up,
    #' write code to wrap up fitting information (from par)
    
    
    pen_deviance_grad_man <- numDeriv::grad(find_pen_deviance, par, sp = sp, y = data$y,
                                            row_list = row_list, basis = basis, k = k)
    pen_deviance_grad_auto <- attr(pen_deviance, "gradient")
    expect_equal(pen_deviance_grad_auto, pen_deviance_grad_man)

    pen_deviance_hess_man <- numDeriv::hessian(find_pen_deviance, par, sp = sp, y = data$y,
                                               row_list = row_list, basis = basis, k = k)
    pen_deviance_hess_auto <- attr(pen_deviance, "hessian")
    expect_equal(pen_deviance_hess_auto, pen_deviance_hess_man)
})
