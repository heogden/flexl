test_that("can fit normal model with fixed k and penalty pars", {
    data_full <- generate_test_data_1()
    data1 <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    sp <- 1000
    sigma <- 0.1
    nbasis <- 10

    mod <- fit_given_sp(data1, sp, 3, nbasis)
    mod10 <- fit_given_sp(data1, 10, 3, nbasis) 
    mod1e5 <- fit_given_sp(data1, 1e5, 3, nbasis)

    sapply(mod10, "[[", "log_ml")
    sapply(mod, "[[", "log_ml")
    sapply(mod1e5, "[[", "log_ml")

    #' now always seems to prefer the smaller tuning par!
    
    plot(data1$x, data1$y)
    curve(mod[[1]]$f0(x), col = 2, add = TRUE, lwd = 2)
    curve(mu, add = TRUE, lty = 2)

    #' fit reasonable (overfitting decreases with sp)

    x_grid <- seq(0, 1, length.out = 100)
    
    curve(mod[[4]]$f(x)[,1], from = 0, to = 1)
    lines(x_grid, sqrt(0.5) * delta(x_grid)[,1], lty = 2)

    curve(mod[[4]]$f(x)[,2], from = 0, to = 1)
    lines(x_grid, sqrt(0.1) * delta(x_grid)[,2], lty = 2)

    sigmas <- sapply(mod, "[[", "sigma")
    
    #' could translate to fraction of "non-residual" variance explained?
    resid_var <- min(sigmas^2)
    non_resid_var <- sigmas^2 - resid_var
    1 - non_resid_var / non_resid_var[1]

    #' could use some cut off in FVE (e.g. 0.99) to select k

    library(tidyverse)
    
    pred_data <- bind_cols(x = data1$x, c = data1$c, eta = eta) %>%
        group_by(c) %>%
        mutate(eta_hat_0 = predict_flexl(mod[[1]], newdata = list(x = x, c = c[1])),
               eta_hat_1 = predict_flexl(mod[[2]], newdata = list(x = x, c = c[1])),
               eta_hat_2 = predict_flexl(mod[[3]], newdata = list(x = x, c = c[1])))

    pred_data_long <- pred_data %>%
        pivot_longer(cols = starts_with("eta_hat_"), names_prefix = "eta_hat_", names_to = "k",
                     values_to = "eta_hat") 
    
    error_k <- pred_data_long %>%
        group_by(k) %>%
        mutate(se = (eta_hat - eta)^2) %>%
        summarise(rmse = sqrt(mean(se)))

    expect_lt(error_k$rmse[2], error_k$rmse[1])
    expect_lt(error_k$rmse[3], error_k$rmse[2])
    
    pred_data_long %>%
    filter(c <= 12) %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = eta_hat, colour = k)) +
    geom_line(aes(y = eta), linetype = "dashed") +
    facet_wrap(vars(c))
    
})
