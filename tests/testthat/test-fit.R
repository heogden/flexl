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
    mod <- fit_flexl(data_norm, nbasis = nbasis, tFVE = 0.999)

    x_pred_data <- crossing(x = seq(from = min(data_norm$x),
                                    to = max(data_norm$x),
                                    length.out = 100),
                            c = unique(data_norm$c))

    pred_data <- x_pred_data  %>%
        group_by(c) %>%
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c[1])))
    
    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat)) +
        geom_point(aes(x = x, y = y), data = data_norm) +
        facet_wrap(vars(c))

    #' test out manual changes to sp
    basis <- find_orthogonal_spline_basis(nbasis, data_norm$x)
    mod_sp <- fit_given_sp_init(data_norm, 100, 10, basis, fve_threshold = 0.99)
    mod_sp$k
    pred_data <- x_pred_data  %>%
        group_by(c) %>%
        mutate(mu_hat = predict_flexl(mod_sp, newdata = list(x = x, c = c[1])))
    
    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat)) +
        geom_point(aes(x = x, y = y), data = data_norm) +
        facet_wrap(vars(c))

    mod_sp$log_ml
    mod$log_ml
    #' k = 4, sp = 100 has slightly higher marginal likelihood
    #' than k = 3, sp = 786

    mod_k4 <- fit_flexl(data_norm, nbasis = nbasis, tFVE = 1, kmax = 4)
    #' why does this fail to give a higher value of log_ml?


    log_ml_poss_a <- c(-266.4289, -234.6486, -204.9176, -181.9478, -171.4674, -160.2171, -158.8954, -158.3792, -168.8195)
    

    mod_k4_sp_man <- fit_given_sp_init(data_norm, mod_k4$sp, 4, basis, fve_threshold = 1)
    
    mod_k4$log_ml
    mod_k4_sp_man$log_ml
    #' there is a difference here

    mod_k4$l_pen
    mod_k4_sp_man$l_pen
    #' two are finding the same optimum value

    plot(mod_k4$par, mod_k4_sp_man$par)
    #' there are some differences in optimum parameters

    #' is there a difference in the fits themselves?
    pred_data <- x_pred_data  %>%
        group_by(c) %>%
        mutate(mu_hat = predict_flexl(mod_k4, newdata = list(x = x, c = c[1])),
               mu_hat_man = predict_flexl(mod_k4_sp_man, newdata = list(x = x, c = c[1])))

     pred_data %>%
        ggplot(aes(x = x)) +
         geom_line(aes(y = mu_hat)) +
         geom_line(aes(y = mu_hat_man), linetype = "dashed") +
         geom_point(aes(x = x, y = y), data = data_norm) +
         facet_wrap(vars(c))
    #' the actual fits look the same
    


})
