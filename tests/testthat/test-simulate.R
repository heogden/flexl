test_that("simulation from fitted model works", {
    data_full <- generate_test_data_1()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    mod <- fit_flexl(data)

    set.seed(1)
    data_sim <- simulate_flexl(mod)

    mod_sim <- fit_flexl(data_sim)
    expect_equal(mod_sim$k, 2)

    expect_gt(cor(mod_sim$par, mod$par), 0.5)

})
