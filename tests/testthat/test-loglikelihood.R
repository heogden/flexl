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

    loglikelihood_grad_man <- numDeriv::grad(find_loglikelihood_k, alpha_k,
                                             X_k = X_k, fit_km1 = fit_km1)
    loglikelihood_grad_auto <- loglikelihood_k_grad(alpha_k, X_k, fit_km1)

    expect_equal(loglikelihood_grad_auto, loglikelihood_grad_man)

    loglikelihood_hess_man <-  numDeriv::hessian(find_loglikelihood_k, alpha_k,
                                                 X_k = X_k, fit_km1 = fit_km1)
    loglikelihood_hess_auto <- loglikelihood_k_hess(alpha_k, X_k, fit_km1)

    expect_equal(loglikelihood_hess_auto, loglikelihood_hess_man)
})
