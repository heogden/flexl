test_that("log_ml contribs are correct", {
    n <- 10
    
    Sigma <- matrix(0.2, nrow = n, ncol = n) + diag(n)
    d <- nrow(Sigma)
    
    log_ml_contrib_man <- -mvtnorm::dmvnorm(rep(0, d), sigma = Sigma, log = TRUE)

    Sigma_inv <- solve(Sigma)
    log_ml_contrib <- approx_log_ml_contrib(Sigma_inv)

    expect_equal(log_ml_contrib, log_ml_contrib_man)
})
