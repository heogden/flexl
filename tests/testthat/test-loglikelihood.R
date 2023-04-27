test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    sp <- 100
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    
    fits <- list()
    fits[[1]] <- fit_0(data, sp, basis)
    fits[[2]] <-  fit_given_k(data, sp, 1, fits[[1]], basis)

    k <- 2
    fit_km1 <- fits[[2]]
    transform <- find_orthogonal_complement_transform(fit_km1$beta)
    X_k <- basis$X %*% transform
    S_k <- emulator::quad.form(basis$S, transform)

    set.seed(1)
    alpha_k <- rnorm(ncol(X_k))

    pen_deviance <- find_pen_deviance_k(alpha_k, sp, X_k, S_k, fit_km1)
    
    pen_deviance_grad_man <- numDeriv::grad(find_pen_deviance_k, alpha_k,
                                            sp = sp, X_k = X_k, S_k = S_k, fit_km1 = fit_km1)
    pen_deviance_grad_auto <- attr(pen_deviance, "gradient")
    expect_equal(pen_deviance_grad_auto, pen_deviance_grad_man)

    pen_deviance_hess_man <- numDeriv::hessian(find_pen_deviance_k, alpha_k,
                                               sp = sp, X_k = X_k, S_k = S_k, fit_km1 = fit_km1)
    pen_deviance_hess_auto <- attr(pen_deviance, "hessian")
    expect_equal(pen_deviance_hess_auto, pen_deviance_hess_man)
})
