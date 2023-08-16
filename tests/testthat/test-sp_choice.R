test_that("log_ml contribs are correct", {
    n <- 10
    
    Sigma <- matrix(0.2, nrow = n, ncol = n) + diag(n)
    d <- nrow(Sigma)
    
    log_ml_contrib_man <- -mvtnorm::dmvnorm(rep(0, d), sigma = Sigma, log = TRUE)

    Sigma_inv <- solve(Sigma)
    log_ml_contrib <- approx_log_ml_contrib(Sigma_inv)

    expect_equal(log_ml_contrib, log_ml_contrib_man)
})

test_that("matches log ML formula for gam in k = 0 case", {
    data_full <- generate_test_data_2()
    data <- data_full$data

    sp <- 100
    nbasis <- 10

    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    S <- basis$S

    M <- 2


    #' replicate (6.21) from Wood (2017), with error fixed
    find_log_ml_0 <- function(sp) {
        mod <- fit_given_sp_init(data, sp, 0, basis)
        sigma2 <- exp(2 * mod$lsigma)
        phi <- sigma2
        S_lop <- S * sp / phi
    
        H <- mod$hessian
        beta_hat <- mod$beta0

        l_p_hat <- mod$l
    
        mod$l_pen + log_det_gen(S_lop, nbasis - 2) / 2 - log_det(H) / 2 + M/2 * log(2 * pi)
    }


    sp_poss <- seq(0.1, 20, length.out = 100)
    log_ml_0_poss <- sapply(sp_poss, find_log_ml_0)
    
    plot(sp_poss, log_ml_0_poss, type = "l")
    #' seems reasonable
    

})
