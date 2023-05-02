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

    sp <- 1000
    sigma <- 0.1
    nbasis <- 10
    
    mod <- fit_given_sp(data1, sp, 3, nbasis)
    
    library(mgcv)
    mod_gam <- gam(y ~ s(x), data = data1)

    sp_poss <- exp(seq(-1, 10, length.out = 20))

    fit_mod <- function(sp) {
        mod <- fit_given_sp(data1, sp, 0, nbasis)
        mod[[1]]
    }
    
    
    mod_poss <- lapply(sp_poss, fit_mod)
    ml_poss <- sapply(mod_poss, "[[", "log_ml")
    
    plot(log(sp_poss), ml_poss, type = "b")
    #' stabilises around 2000 or so 
    sp_poss[which.max(ml_poss)]
    
    sapply(mod10, "[[", "log_ml")
    sapply(mod, "[[", "log_ml")
    sapply(mod1e5, "[[", "log_ml")

    mod10[[1]]$lpen_hat
    mod1e5[[1]]$lpen_hat
    
    #' now always seems to prefer the smaller tuning par!

    data <- data1
    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    Matrix::rankMatrix(basis$S)
    #' a rank 8 matrix (while dimension is 10)
    
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
