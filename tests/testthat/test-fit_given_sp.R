test_that("can fit normal model with fixed k and penalty pars", {
    data_full <- generate_test_data_1()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    sp <- 100
    
    system.time(mod <- fit_given_sp(data, sp, 3, nbasis, 0.99))

    library(tidyverse)
    
    pred_data <- bind_cols(x = data$x, c = data$c, eta = eta) %>%
        group_by(c) %>%
        mutate(eta_hat_0 = predict_flexl(mod[[1]], newdata = list(x = x, c = c[1])),
               eta_hat_1 = predict_flexl(mod[[2]], newdata = list(x = x, c = c[1])),
               eta_hat_2 = predict_flexl(mod[[3]], newdata = list(x = x, c = c[1])),
               eta_hat_3 = predict_flexl(mod[[4]], newdata = list(x = x, c = c[1])))

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

