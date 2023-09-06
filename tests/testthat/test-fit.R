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
        mutate(eta_hat = predict_flexl(mod, newdata = list(x = x, c = c)))

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

    #' check prediction
    y_hat_data <- predict_flexl(mod, newdata = data)
    fitted_y <- fitted_flexl(mod)
    expect_equal(y_hat_data, fitted_y)

    
    mu_2_fun <- function(x) {
        predict_flexl(mod, newdata = data.frame(x = x, c = 2))
    }

    newdata <- data.frame(x = seq(min(data$x), max(data$x), length = 10),
                          c = 2)
    d_mu_hat_data_man <- numDeriv::grad(mu_2_fun, newdata$x)
    d_mu_hat_data <- predict_flexl(mod, newdata = newdata, deriv = TRUE)
    expect_equal(d_mu_hat_data, d_mu_hat_data_man)
    
    #' look at uncertainty
    n_samples <- 1000
    samples <- find_samples(mod, n_samples)

    x_pred_data <- crossing(x = seq(from = min(data$x),
                                    to = max(data$x),
                                    length.out = 100),
                            c = unique(data$c))

    pred_data <- x_pred_data  %>%
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c), interval = TRUE,
                                      samples = samples),
               d_mu_hat = predict_flexl(mod, newdata = list(x = x, c = c), deriv = TRUE,
                                        interval = TRUE, samples = samples)) %>%
        group_by(c) %>%
        mutate(mu = data_full$eta_fun(x, c[1]),
               d_mu = numDeriv::grad(data_full$eta_fun, x = x, c = c[1]))

    pred_data %>%
        filter(c <= 12) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat$estimate)) +
        geom_line(aes(y = mu), colour = "red", linetype = "dashed") + 
        geom_ribbon(aes(ymin = mu_hat$lower, ymax = mu_hat$upper), alpha = 0.2) + 
        geom_point(aes(x = x, y = y), data = data %>% filter(c <= 12)) +
        facet_wrap(vars(c))


    coverage <- as.numeric(pred_data %>%
        mutate(covers = ((mu_hat$lower < mu) & (mu_hat$upper > mu)),) %>%
        ungroup() %>%
        summarise(coverage = mean(covers)))

    expect_gt(coverage, 0.9)
    expect_lt(coverage, 1)

    
    pred_data %>%
        filter(c <= 12) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = d_mu_hat$estimate)) +
        geom_line(aes(y = d_mu), colour = "red", linetype = "dashed") + 
        geom_ribbon(aes(ymin = d_mu_hat$lower, ymax = d_mu_hat$upper), alpha = 0.2) + 
        facet_wrap(vars(c))

    coverage_d <- as.numeric(pred_data %>%
        mutate(d_covers = ((d_mu_hat$lower < d_mu) & (d_mu_hat$upper > d_mu)),) %>%
        ungroup() %>%
        summarise(d_coverage = mean(d_covers)))

    expect_gt(coverage_d, 0.9)
    expect_lt(coverage_d, 1)

    
    
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
        mutate(eta_hat = predict_flexl(mod, newdata = list(x = x, c = c)))

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
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c)))

    
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

    nbasis <- 15
    mod <- fit_flexl(data, nbasis = nbasis)


    x_pred_data <- crossing(x = seq(from = min(data$x),
                                    to = max(data$x),
                                    length.out = 100),
                            c = unique(data$c))

    pred_data <- x_pred_data  %>%
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c)))
    
    pred_data %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat)) +
        geom_point(aes(x = x, y = y), data = data) +
        facet_wrap(vars(c))
    
    
})


test_that("fits the fat data", {
    library(ALA)
    library(tidyverse)
    
    data <- fat %>%
        as_tibble %>%
        mutate(c = as.integer(as.factor(id)),
               x = time.menarche,
               y = percent.fat,
               .keep = "none")

    mod <- fit_flexl(data)
  
    n_samples <- 1000
    samples <- find_samples(mod, n_samples)

    x_pred_data <- crossing(x = seq(from = min(data$x),
                                    to = max(data$x),
                                    length.out = 100),
                            c = unique(data$c))

    pred_data <- x_pred_data  %>%
        mutate(mu_hat = predict_flexl(mod, newdata = list(x = x, c = c), interval = TRUE,
                                      samples = samples))

    pred_data %>%
        filter(c <= 12) %>%
        ggplot(aes(x = x)) +
        geom_line(aes(y = mu_hat$estimate)) +
        geom_ribbon(aes(ymin = mu_hat$lower, ymax = mu_hat$upper), alpha = 0.2) + 
        geom_point(aes(x = x, y = y), data = data %>% filter(c <= 12)) +
        facet_wrap(vars(c))
    
    set.seed(1)

    data_noise <- data %>%
        mutate(y = y + rnorm(length(y), 0, 20))

    mod_noise <- fit_flexl(data_noise, trace = TRUE)

    expect_lte(mod_noise$k, mod$k)
    expect_gte(mod_noise$sp, mod$sp)
})

test_that("fits the cd4 data", {
    library(refund)
    library(tidyverse)
    
     #' modified from refund::ccb.fpc
    data(cd4)
    #' obtain a subsample of the data with 25 subjects
    set.seed(1236)
    sample = sample(1:dim(cd4)[1], 25)
    Y.sub = cd4[sample,]

    times <- as.numeric(colnames(Y.sub))
    data <- tibble(c = row(Y.sub)[!is.na(Y.sub)],
                   y = Y.sub[!is.na(Y.sub)],
                   x = times[col(Y.sub)[!is.na(Y.sub)]])
    
    y_hat <- fitted_flexl(mod)
    y_hat_8 <- y_hat[data$c == 8]
    y_8 <- data$y[data$c == 8]

    #' gives very bad fit for cluster 8 at the moment
    expect_false(all(y_hat_8 > y_8))

  
})
