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

    expect_gt(mod$sp, 0.1)
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

    expect_gt(mod$sp, 0.1)


    nbasis <- 30
    tFVE <- 0.99
    
    lsp_poss <- seq(-5, 10, length.out = 10)
    sp_poss <- exp(lsp_poss)
    
    nbasis <- 30
    kmax <- 2
    sp <- sp_poss[1]
    
    fit_sp1 <- fit_given_sp(data, sp, kmax, nbasis)

    fit_sp1[[1]]$opt
    ## Problem here: didn't previously converge, have now upped maxit
    fit_sp1[[2]]$opt
    ## Continue from here: do we now get better log_ml 

    

    

    fits_poss <- list()
    fits_poss[[1]] <- fit_given_sp(data, sp_poss[1], 10, nbasis, tFVE)
    kmax <- length(fits_poss[[1]]) - 1

    log_ml_poss <- c()
    log_ml_poss[1] <- (fits_poss[[1]])[[kmax + 1]]$log_ml
    
    for(i in 2:length(sp_poss)) {
        fits_poss[[i]] <- fit_given_sp(data, sp_poss[i], kmax, nbasis, 1)
        log_ml_poss[i] <- (fits_poss[[i]])[[kmax + 1]]$log_ml
        #f(i > 3 & (log_ml_poss[i] < log_ml_poss[i-1])) {
        #   lsp_poss <- lsp_poss[1:i]
        #   break
        #
    }

    plot(lsp_poss, log_ml_poss, type = "l")
    #' what happens if we try to use parameter estimates from
    #' first sp as starting points?

    sp <- sp_poss[2]

    fit_orig <- fits_poss[[2]][[length(fits_poss[[2]])]]

    fit_sp1 <-  fits_poss[[1]][[length(fits_poss[[1]])]]
    basis <- find_orthogonal_spline_basis(nbasis, data$x)
    
    system.time(fit_without_sp1 <- fit_given_sp(data, sp, kmax, nbasis, 1))
    

    system.time(fit_using_sp1 <- fit_given_fit_other_sp(data, sp, fit_sp1, basis))

    fit_sp4 <-  fits_poss[[4]][[length(fits_poss[[4]])]]
    system.time(fit_using_sp4 <- fit_given_fit_other_sp(data, sp_poss[4] + 1, fit_sp4, basis))

    
    fit_without_sp1[[3]]$opt
    fit_using_sp1$opt
    #' still needs to do a lot of runs, even using sp1
    #' would using hessian information help?
    #' could get better starting point?
    g_sp2 <-  loglikelihood_pen_grad(fit_sp1$par, basis$X, data$y, data$c - 1, sp, basis$S, kmax)

    par_diff_pred <- -solve(H, g_sp2)

    par_diff_real <- fit_using_sp1$par - fit_sp1$par

    plot(par_diff_pred, par_diff_real, xlim = c(-0.1, 0.02), ylim = c(-0.1, 0.02))
    abline(a = 0, b = 1)

    par_pred <- fit_sp1$par + par_diff_pred
    fit_par_pred <- fit_sp1
    fit_par_pred$par <- par_pred
    
    system.time(fit_using_par_pred <- fit_given_fit_other_sp(data, sp, fit_par_pred, basis))
    fit_using_par_pred$opt

    system.time(H <- loglikelihood_pen_hess(fit_sp1$par, basis$X, data$y, data$c - 1, sp, basis$S, kmax))
    

    
    
    fit_orig$l_pen
    fit_using_sp1$l_pen

    #' they have different penalised loglikelihoods.
    #' (but now quite similar)
    #' probably reduced costs a lot in using other fit: check this

    curve(fit_using_sp1$f0(x), from = 0, to = 1)
    curve(fit_orig$f0(x), lty = 2, add = TRUE)
    #' something could be going wrong with estimating mean function?
    #' actually here, true mean is zero, so not clear which one is wrong

    curve(fit_using_sp1$f(x)[,1], from = 0, to = 1)
    curve(f1, lty = 2, add = TRUE)
    #' we don't have scaling here

    curve(fit_orig$f(x)[,1], from = 0, to = 1)
    curve(f1, lty = 2, add = TRUE)
    #' we don't have right scaling here?


    curve(fit_orig$f(x)[,1], from = 0, to = 1)
    curve(fit_using_sp1$f(x)[,1], from = 0, to = 1, col = 2, add = TRUE)
    #' very different scales.

    sd(fit_orig$u[,1])
    sd(fit_orig$u[,2])
    cor(fit_orig$u[,1], fit_orig$u[,2])
    
    sd(fit_using_sp1$u[,1])
    sd(fit_using_sp1$u[,2])
    cor(fit_using_sp1$u[,1], fit_using_sp1$u[,2])

    #' These should be close to 1, 1, 0: something has gone wrong here.

    
    
    
})
