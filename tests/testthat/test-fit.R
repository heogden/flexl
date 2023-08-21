test_that("sensible fit for test data 1 (straight lines)", {
    data_full <- generate_test_data_1()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta

    mod <- fit_flexl(data)

    expect_equal(mod$k, 2)
    expect_gt(mod$sp, 1000)

    library(tidyverse)
    
    pred_data <- bind_cols(x = data$x, c = data$c, eta = eta) %>%
        group_by(c) %>%
        mutate(eta_hat = predict_flexl(mod, newdata = list(x = x, c = c[1])))

    rmse <- sqrt(mean(pred_data$eta_hat - pred_data$eta)^2)
    expect_lt(rmse, 0.1)

    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = eta_hat)) +
        geom_line(aes(y = eta), linetype = "dashed") +
        facet_wrap(vars(c))

    sd(mod$u[,1])
    sd(mod$u[,2])
    cor(mod$u[,1], mod$u[,2])
})


test_that("sensible fit for test data 2 (not straight lines)", {
    data_full <- generate_test_data_2()
    data <- data_full$data
    mu <- data_full$mu
    delta <- data_full$delta
    eta <- data_full$eta


    data %>%
        filter(c <= 10) %>%
        ggplot(aes(x = x, y = y)) +
        geom_point() +
        facet_wrap(vars(c))
    
    mod <- fit_flexl(data)

    expect_equal(mod$k, 2)
    expect_lt(mod$sp, 1000)

    library(tidyverse)
    
    pred_data <- bind_cols(x = data$x, c = data$c, eta = eta) %>%
        group_by(c) %>%
        mutate(eta_hat = predict_flexl(mod, newdata = list(x = x, c = c[1])))

    rmse <- sqrt(mean(pred_data$eta_hat - pred_data$eta)^2)
    expect_lt(rmse, 0.1)

    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = eta_hat)) +
        geom_line(aes(y = eta), linetype = "dashed") +
        facet_wrap(vars(c))
})

test_that("chooses a reasonably large sp if have straight lines", {
    sigma_u <- 0.5
    sigma <- 0.1
    n_clusters <- 50
    n_obs_per_cluster <- 10
    beta_0 <- -1
    beta_1 <- 2

    set.seed(1)

    u <- rnorm(n_clusters, sd = sigma_u)

    x_full <- seq(0, 1, length.out = 1000)
    x_index <- sample(seq_along(x_full), n_clusters * n_obs_per_cluster, replace = TRUE)

    library(tidyverse)

    data <- bind_cols(x_index = x_index,
                      c = rep(1:n_clusters, each = n_obs_per_cluster)) %>%
        mutate(x = x_full[x_index],
               eta = (beta_0 + beta_1 * x + u[c]),
               epsilon = rnorm(length(x_index), 0, sigma),
               y = eta + epsilon)

    eta_data <- crossing(x = x_full,
                         c = 1:n_clusters) %>%
        mutate(eta = (beta_0 + beta_1 * x + u[c]))

    mod <- fit_flexl(data)

    expect_gt(mod$sp, 100)
})



test_that("error if data has missing values", {
    library(ALA)
    library(tidyverse)
    
    data <-  exercise %>%
        as_tibble %>%
        mutate(c = as.integer(as.factor(id)),
               x = day,
               y = strength,
               .keep = "none")

    expect_error(fit_flexl(data), "missing values")

})


test_that("gives reasonable fit with tricky blip function", {
    
    g1 <- function(x) {
        dnorm(x, mean = 0.5, sd = 0.05)
    }


    const <- integrate(function(x){g1(x)^2}, lower = 0, upper = 1)$value

    f1 <- function(x) {
        g1(x) / sqrt(const)
    }
    #' so f1 is a normalised version of g


    #' Next, generate data
    
    set.seed(1)

    n <- 100
    u <- rnorm(n)
    n_i <- 10
    sigma <- 0.01
    
    subject <- rep(1:n, each = n_i)
    t <- runif(n * n_i, 0, 1)

    mu <- u[subject] * f1(t)

    epsilon <- rnorm(n * n_i, sd = sigma)
    y <- mu + epsilon

    data <- data.frame(x = t,
                       y = y,
                       c = subject)

    
    mod <- fit_flexl(data, nbasis = 30)


    library(tidyverse)

    pred_data<- crossing(x = seq(0, 1, length.out = 100),
                         c = 1:n) %>%
        mutate(mu = u[c] * f1(x)) %>%
        group_by(c) %>%
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c[1])))

    
    rmse <- sqrt(mean(pred_data$mu_hat - pred_data$mu)^2)
    expect_lt(rmse, 0.1)

    pred_data %>%
        filter(c <= 20) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat)) +
        geom_line(aes(y = mu), linetype = "dashed") +
        facet_wrap(vars(c))
    
    
})

test_that("fits the sleepstudy data", {
    library(lme4)
    library(tidyverse)
    
    data <- sleepstudy %>%
        as_tibble %>%
        mutate(c = as.integer(as.factor(Subject)),
               y = Reaction,
               x = Days)

    data_norm <- data
    data_norm$y <- (data$y - mean(data$y)) / sd(data$y)
    data_norm$x <- (data$x - mean(data$x)) / sd(data$x)

    #' TODO: automatically normalise data for numerical stability.
    nbasis <- 15
    mod <- fit_flexl(data_norm, nbasis = nbasis)
    
    mod_2 <- fit_flexl(data_norm, nbasis = nbasis, lsp_poss = 5:15)
    #' Two lsp_poss now choose same sp, k
    #' BUT different fits
    mod$l_pen
    split_alpha(mod$alpha, nbasis, mod$k) #' all of the fourth component v small
    mod$log_ml
    correct_lprior_alpha(mod, basis)
    log_det(mod$hessian)
    
    mod_2$l_pen
    split_alpha(mod_2$alpha, nbasis, mod_2$k) #' other components are also quite different
    #' (but fit looks the same in the two cases overall)
    
    mod_2$log_ml
    correct_lprior_alpha(mod_2, basis)
    log_det(mod_2$hessian)

    x_pred_data <- crossing(x = seq(from = min(data_norm$x),
                                    to = max(data_norm$x),
                                    length.out = 100),
                            c = unique(data_norm$c))

    pred_data <- x_pred_data  %>%
        group_by(c) %>%
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c[1])),
               mu_hat_2 = predict_flexl(mod_2, newdata = list(x = x, c = c[1])))
    )
    
    pred_data %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = mu_hat)) +
    geom_line(aes(y = mu_hat_2), lty = 2) +
    geom_point(aes(x = x, y = y), data = data_norm) +
    facet_wrap(vars(c))
    
    sp <- exp(-5)

    
    basis <- find_orthogonal_spline_basis(nbasis, data_norm$x)


    find_log_ml_k_given_sp <- function(sp, kmax) {
        log_ml_poss <- c()
        mod <- fit_0(data_norm, sp, basis)
        log_ml_poss[1] <- mod$log_ml
        for(k in 1:kmax) {
            mod_prev <- mod
            mod <- fit_given_fit_km1(data_norm, sp, k, mod_prev, basis)
            log_ml_poss[k+1] <- mod$log_ml
        }
        log_ml_poss
    }

    kmax <- 6
    log_ml_sp_1 <- find_log_ml_k_given_sp(exp(-5), kmax)
    log_ml_sp_2 <- find_log_ml_k_given_sp(exp(-4), kmax)
    log_ml_sp_3 <- find_log_ml_k_given_sp(exp(-3), kmax)
    log_ml_sp_4 <- find_log_ml_k_given_sp(exp(-2), kmax)
    log_ml_sp_5 <- find_log_ml_k_given_sp(exp(-1), kmax)
    log_ml_sp_6 <- find_log_ml_k_given_sp(exp(0), kmax)
    log_ml_sp_7 <- find_log_ml_k_given_sp(exp(1), kmax)

    plot(0:kmax, log_ml_sp_2, type = "b")
    plot(0:kmax, log_ml_sp_4, type = "b")
    plot(0:kmax, log_ml_sp_5, type = "b")
    plot(0:kmax, log_ml_sp_6, type = "b")
    plot(0:kmax, log_ml_sp_7, type = "b")


    #' what happens if we compare by using FVE instead?

    kmax <- 10
    
    sp <- exp(-5)
    mods_m5 <- list(fit_0(data_norm, sp, basis))
    for(k in 1:kmax) {
        mods_m5[[k+1]] <- fit_given_fit_km1(data_norm, sp, k, mods_m5[[k]], basis)
    }

    sapply(mods_m5, "[[", "sigma")
    #' sigma still substantially reducing, all the way to k = 9. Better to use FVE instead.

    sp <- 1
    mods_0 <- list(fit_0(data_norm, sp, basis))
    for(k in 1:kmax) {
        mods_0[[k+1]] <- fit_given_fit_km1(data_norm, sp, k, mods_0[[k]], basis)
    }

    sapply(mods_0, "[[", "sigma")
    #' gives k = 7

    sp <- exp(6)
    mods_6 <- list(fit_0(data_norm, sp, basis))
    for(k in 1:kmax) {
        mods_6[[k+1]] <- fit_given_fit_km1(data_norm, sp, k, mods_6[[k]], basis)
    }

    sigmas <- sapply(mods_6, "[[", "sigma")
    log_mls <- sapply(mods_6, "[[", "log_ml")
    r <- sigmas[-1]/sigmas[-length(sigmas)]
    which(r > 1- 1e-3)
    #' can drop k = 4 and above

    #' same conclusions looking at lambdas
    colMeans((mod$f(data_norm$x)^2))
    #' contrib from k = 4 is very small
    

    #' This gives us a way to choose k, for fixed sigma
    #' But how should we compare sigma?

    mod_m5 <- mods_m5[[10]]
    mod_0 <- mods_0[[8]]
    mod_6 <- mods_6[[5]]

    
})

