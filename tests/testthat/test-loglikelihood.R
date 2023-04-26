test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()
    
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

    pen_loglikelihood <- find_pen_loglikelihood_k(alpha_k, sp, X_k, S_k, fit_km1)
    
    pen_loglikelihood_grad_man <- numDeriv::grad(find_pen_loglikelihood_k, alpha_k,
                                                 sp = sp, X_k = X_k, S_k = S_k, fit_km1 = fit_km1)
    pen_loglikelihood_grad_auto <- attr(pen_loglikelihood, "gradient")
    expect_equal(pen_loglikelihood_grad_auto, pen_loglikelihood_grad_man)

    pen_loglikelihood_hess_man <- numDeriv::hessian(find_pen_loglikelihood_k, alpha_k,
                                                    sp = sp, X_k = X_k, S_k = S_k, fit_km1 = fit_km1)
    pen_loglikelihood_hess_auto <- attr(pen_loglikelihood, "hessian")
    expect_equal(pen_loglikelihood_hess_auto, pen_loglikelihood_hess_man)
})
