test_that("fitting with mean-only works correctly", {
    data_full <- generate_test_data_1()
    data <- data_full$data

    nbasis <- 10
    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    sp <- 100
    mod <- fit_0(data, sp, basis)

    X <- basis$X
    S <- basis$S

    l_pen_beta <- function(beta, sigma) {
        l <- sum(dnorm(data$y, X %*% beta, sd = sigma, log = TRUE))
        spr <- sp / (2 * sigma^2)
        pen <- - spr * emulator::quad.form(S, beta)
        l + pen
    }

    opt_man <- optim(mod$beta_0, l_pen_beta, sigma = mod$sigma, control = list(fnscale = -1),
                     method = "BFGS")

    expect_equal(mod$beta_0, opt_man$par)
    
})


test_that("can fit normal model with fixed k and penalty pars", {
    data_full <- generate_test_data_1()
    data1 <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    nbasis <- 10

    t1000 <- system.time(mod <- fit_given_sp(data1, 1000, 3, nbasis))
    t100 <-  system.time(mod100 <- fit_given_sp(data1, 100, 3, nbasis))
    t100_other_fit <- system.time(mod100_other_fit <- fit_given_sp(data1, 100, 3, nbasis, fit_other_sp = mod))

    expect_lt(t100_other_fit[1], t100[1])
    
    expect_equal(mod100_other_fit[[4]]$beta, mod100[[4]]$beta, tolerance = 1e-4)
    expect_equal(mod100_other_fit[[4]]$l_pen_hat, mod100[[4]]$l_pen_hat)

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

test_that("can fit with full-loglikelihood stage", {
    data_full <- generate_test_data_1()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    nbasis <- 10
    k <- 2
    sp <- 100
    
    mod <- fit_given_sp(data, sp, k, nbasis, full = TRUE)
    mod_not_full <- fit_given_sp(data, sp, k, nbasis, full = FALSE)

    pred_data <- bind_cols(x = data$x, c = data$c, eta = eta) %>%
        group_by(c) %>%
        mutate(eta_hat_full = predict_flexl(mod[[k+1]], newdata = list(x = x, c = c[1])),
               eta_hat_not_full = predict_flexl(mod_not_full[[k+1]], newdata = list(x = x, c = c[1])))

    pred_data_long <- pred_data %>%
        pivot_longer(cols = starts_with("eta_hat_"), names_prefix = "eta_hat_", names_to = "type",
                     values_to = "eta_hat") 
    
    error_type <- pred_data_long %>%
        group_by(type) %>%
        mutate(se = (eta_hat - eta)^2) %>%
        summarise(rmse = sqrt(mean(se)))

    expect_lte(error_type$rmse[1], error_type$rmse[2])

    pred_data_long %>%
        filter(c <= 12) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = eta_hat, colour = type)) +
        geom_line(aes(y = eta), linetype = "dashed") +
        geom_point(aes(x = x, y = y), data = (data %>% filter(c <= 12))) +
        facet_wrap(vars(c))
    
})
