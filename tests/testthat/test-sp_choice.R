test_that("log_ml contribs are correct", {
    n <- 10
    
    Sigma <- matrix(0.2, nrow = n, ncol = n) + diag(n)
    d <- nrow(Sigma)
    
    log_ml_contrib_man <- -mvtnorm::dmvnorm(rep(0, d), sigma = Sigma, log = TRUE)

    Sigma_inv <- solve(Sigma)
    log_ml_contrib <- approx_log_ml_contrib(Sigma_inv)

    expect_equal(log_ml_contrib, log_ml_contrib_man)
})

test_that("Log ML not always increasing with k", {
    
    find_log_ml_curr <- function(sp, k, data, nbasis) {
        basis <- find_orthogonal_spline_basis(nbasis, data$x)
        S <- basis$S
        
        mod <- fit_given_sp_init(data, sp, k, basis)
        mod$log_ml
    }

    find_log_ml <- function(sp, k, data, nbasis) {
        basis <- find_orthogonal_spline_basis(nbasis, data$x)
        S <- basis$S
        
        mod <- fit_given_sp_init(data, sp, k, basis)
        
        sigma2 <- exp(2 * mod$lsigma)
        phi <- sigma2
        S_lop <- S * sp / phi
    
        H <- mod$hessian
        beta_hat <- mod$beta0

        l_p_hat <- mod$l

        r <- nbasis - 2
        
        cat("terms are: \n")
        cat("l_pen = ", mod$l_pen, "\n")
        cat("correction = ", (k+1) * log_det_gen(S_lop, r) / 2 - (k+1)*r/2 * log(2 * pi), "\n")
        cat("p/2 * log(2*pi) = ",  (length(mod$par))/2 * log(2 * pi), "\n")
        cat("1/2 * log_det(-hessian) = ", - log_det(H)/2, "\n")
        
        mod$l_pen - log_det(H)/2 +
            (k+1) * log_det_gen(S_lop, r) / 2  +
            (length(mod$par) - (k+1)*r)/2 * log(2 * pi)
    }

    nbasis <- 10
    
    #' generate data with truth k = 0
    data0 <- generate_test_data_0()$data

    sp_poss <- seq(0.01, 1, length.out = 100)
    log_ml_0_poss <- sapply(sp_poss, find_log_ml, k = 0, data = data0, nbasis = nbasis)
    log_ml_1_poss <- sapply(sp_poss, find_log_ml, k = 1, data = data0, nbasis = nbasis)

    plot(range(sp_poss), range(c(log_ml_0_poss, log_ml_1_poss)), type = "n")
    lines(sp_poss, log_ml_0_poss)
    lines(sp_poss, log_ml_1_poss, lty = 2)

    expect_gt(max(log_ml_0_poss), max(log_ml_1_poss))

    log_ml_0_poss_curr <- sapply(sp_poss, find_log_ml_curr, k = 0, data = data0, nbasis = nbasis)
    log_ml_1_poss_curr <- sapply(sp_poss, find_log_ml_curr, k = 1, data = data0, nbasis = nbasis)

    plot(range(sp_poss), range(c(log_ml_0_poss_curr, log_ml_1_poss_curr)), type = "n")
    lines(sp_poss, log_ml_0_poss_curr)
    lines(sp_poss, log_ml_1_poss_curr, lty = 2)

    expect_gt(max(log_ml_0_poss_curr), max(log_ml_1_poss_curr))

    plot(range(sp_poss), range(c(log_ml_0_poss, log_ml_0_poss_curr)), type = "n")
    lines(sp_poss, log_ml_0_poss)
    lines(sp_poss, log_ml_0_poss_curr, lty = 2)
    #' two results the same!


    plot(range(sp_poss), range(c(log_ml_1_poss, log_ml_1_poss_curr)), type = "n")
    lines(sp_poss, log_ml_1_poss)
    lines(sp_poss, log_ml_1_poss_curr, lty = 2)
    #' two results different, but similar. Where is the difference?


    #' Look at one sp
    sp <- 2
    find_log_ml_curr(sp, 0, data0, nbasis)
    find_log_ml(sp, 0, data0, nbasis)
    #' next: look at individual terms here
    
    #' generate data with truth k = 2
    data2 <- generate_test_data_2()$data
    
    sp_poss <- seq(0.01, 10, length.out = 100)
    log_ml_0_poss <- sapply(sp_poss, find_log_ml, k = 0, data = data2, nbasis = nbasis)
    log_ml_1_poss <- sapply(sp_poss, find_log_ml, k = 1, data = data2, nbasis = nbasis)

    plot(range(sp_poss), range(c(log_ml_0_poss, log_ml_1_poss)), type = "n")
    lines(sp_poss, log_ml_0_poss)
    lines(sp_poss, log_ml_1_poss, lty = 2)

    expect_gt(max(log_ml_1_poss), max(log_ml_0_poss))


    sp <- 2

    nbasis <- 10
    basis <- find_orthogonal_spline_basis(nbasis, data2$x)
    S <- basis$S
    
    mod <- fit_given_sp_init(data2, sp, 3, basis)
    #' find example T (may not match exactly)

    beta <- mod$beta

    #' finding transform for k
    k <- 3

    T_list <- find_T_list(alpha, nbasis, k)
    T_k <- T_list[[k]]
    
    S_k <- t(T_k) %*% S %*% T_k
    
    alpha <- mod$alpha

    alpha_list <- split_alpha(alpha, nbasis, k)
    alpha_k <- alpha_list[[k]]

    beta_k <- beta[,k]
    qform_beta <- emulator::quad.form(S, beta_k)

    qform_alpha <- emulator::quad.form(S_k, alpha_k)

    expect_equal(qform_alpha, qform_beta)

    S_k_los2 <- sp / mod$sigma^2 * S_k
    log_det_gen(S_k_los2, r_k)
    (log_det_gen(S_k, r_k) + r_k * (log(sp) - 2 * mod$lsigma))
    

    library(Matrix)
    rankMatrix(S)
    rankMatrix(S_k)
    #' rank(S_k) = min(rank(S), nrow(S_k)) = min(n_B - 2, n_B - k + 1)
    #' so rank(S_k) = n_B - 2 for k = 1, 2,
    #'                n_B - (k-1) for k = 3, ....
    #' we need n_B > k-1 for this to make sense
    #' so always take k < n_B + 1, i.e. k <= n_B

    r_k <- min(nbasis - 2, nbasis - k + 1)
    
    log_det_gen(S_k, r_k)
    log_det_gen(S, nbasis - 2)
    #' close if k = 2 (could be numerical error)

    evs_S <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    evs <- eigen(S_k, symmetric = TRUE, only.values = TRUE)$values

    log(evs_S[1:r_k]) - log(evs)
    #' most of the difference comes in last component for k = 3
    
    log_det_gen(S_k, r_k) - log_det_gen(S, nbasis - 2)
})
