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

    ld_auto <- ldmvnorm(a, info_km1)
    ld_man <- mvtnorm::dmvnorm(z, sigma = Sigma_k, log = TRUE)

    expect_equal(ld_auto, ld_man)

    ld_grad_auto <- ldmvnorm_grad(a, info_km1)
    ld_grad_man <- numDeriv::grad(ldmvnorm, a, info_km1 = info_km1)

    expect_equal(ld_grad_auto, ld_grad_man, tolerance = 1e-4)
 

    log_det_Sigma <- function(a, info_km1) {
          b <- info_km1$Sigma_inv %*% a
          c <- 1 + sum(a * b)
          log(c) + info_km1$ldet_Sigma
    }
    
    b <- as.numeric(info_km1$Sigma_inv %*% a)
    
    c <- 1 + sum(a * b)
    dc <- 2 * b
    d2c <- 2 * info_km1$Sigma_inv

    log_det_Sigma_hess_auto <- log_det_Sigma_hess(c, dc, d2c)
    log_det_Sigma_hess_man <- numDeriv::hessian(log_det_Sigma, a, info_km1 = info_km1)

    expect_equal(log_det_Sigma_hess_auto, log_det_Sigma_hess_man, tolerance = 1e-4)

    quad_form <- function(a, info_km1) {
         b <- info_km1$Sigma_inv %*% a
         c <- 1 + sum(a * b)
         ldet_Sigma <- log(c) + info_km1$ldet_Sigma
         d <- sum(info_km1$Sigma_inv_z * a)
         
         info_km1$tz_Sigma_inv_z - d^2/c
    }

        
    d <- sum(info_km1$Sigma_inv_z * a)
    dd <- as.numeric(info_km1$Sigma_inv_z)

    
    quad_form_hess_auto <- quad_form_hess(c, dc, d2c, d, dd)
    quad_form_hess_man <- numDeriv::hessian(quad_form, a, info_km1 = info_km1)
    
    expect_equal(quad_form_hess_auto, quad_form_hess_man, tolerance = 1e-4)

    ld_hess_auto <- ldmvnorm_hess(a, info_km1)
    ld_hess_man <- numDeriv::hessian(ldmvnorm, a, info_km1 = info_km1)

    expect_equal(ld_hess_auto, ld_hess_man, tolerance = 1e-4)
 


})
