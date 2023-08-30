test_that("prediction works with unordered data with repeated x values", {
    library(tidyverse)
    
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

    data <- bind_cols(x_index = x_index,
                      c = rep(1:n_clusters, each = n_obs_per_cluster)) %>%
        mutate(x = x_full[x_index],
               eta = (beta_0 + beta_1 * x + u[c]),
               epsilon = rnorm(length(x_index), 0, sigma),
               y = eta + epsilon)

    mod <- fit_flexl(data)
    y_hat <- fitted_flexl(mod)
    y_hat_pred <- predict_flexl(mod, newdata = data)

    expect_equal(y_hat, y_hat_pred)
})

