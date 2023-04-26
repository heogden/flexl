test_that("ldmvnorm and derivatives work", {

    set.seed(1)
    n <- 100
    
    Sigma_km1 <- matrix(0.2, nrow = n, ncol = n) + diag(n)
    a <- rnorm(n, sd = 0.5)

    Sigma_k <- Sigma_km1 + outer(a, a)
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_k))
    
    Sigma_km1_inv <- solve(Sigma_km1)
    ldet_Sigma_km1 <- determinant(Sigma_km1, logarithm = TRUE)$modulus
    attr(ldet_Sigma_km1, "logarithm") <- NULL

    info_km1 <- list(Sigma_inv = Sigma_km1_inv,
                     ldet_Sigma = ldet_Sigma_km1,
                     z = z,
                     Sigma_inv_z = Sigma_km1_inv %*% z,
                     tz_Sigma_inv_z = emulator::quad.form(Sigma_km1_inv, z))

    ld_auto_full <- ldmvnorm(a, info_km1)

    ld_auto <- drop_attributes(ld_auto_full)
    ld_man <- mvtnorm::dmvnorm(z, sigma = Sigma_k, log = TRUE)

    expect_equal(ld_auto, ld_man)

    ld_grad_auto <- attr(ld_auto_full, "gradient")
    ld_grad_man <- numDeriv::grad(ldmvnorm, a, info_km1 = info_km1)

    expect_equal(ld_grad_auto, ld_grad_man, tolerance = 1e-4)

    ld_hess_auto <- attr(ld_auto_full, "hessian")
    ld_hess_man <- numDeriv::hessian(ldmvnorm, a, info_km1 = info_km1)

    expect_equal(ld_hess_auto, ld_hess_man, tolerance = 1e-4)
})
