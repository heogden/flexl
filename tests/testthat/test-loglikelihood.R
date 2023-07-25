test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    k <- 2
    sp <- 10
    basis <- find_orthogonal_spline_basis(nbasis, data$x)


    row_list <- split(1:nrow(data), data$c)
    alpha_components <- find_alpha_components(nbasis, k)

    set.seed(1)
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), 1)

    lp <- find_loglikelihood_pen_noderiv(par, basis$X, data$y, data$c - 1, sp, basis$S, k)

    l_pen_R <- function(par, sp) {
        par_split <- split_par(par, basis$nbasis)

        beta0 <- par_split$beta0
        beta <- find_beta(par_split$alpha, basis$nbasis, k)$value
        sigma <- exp(par_split$lsigma)
        
        l <- find_loglikelihood_beta(beta0, beta, sigma, data$y, row_list, basis, k)
        beta_full <- cbind(beta0, beta)
        w_each <- apply(beta_full, 2, find_wiggliness_f_k, S = basis$S, derivs = FALSE)
        w <- sum(w_each)
        spr <- sp / (2 * sigma^2)
        pen <- spr * w
        l - pen
    }

    lp_R <- l_pen_R(par, sp)

    library(testthat)
    expect_equal(lp, lp_R)

    lp_with_deriv <- find_loglikelihood_pen_with_hess(par, basis$X, data$y, data$c - 1, sp,
                                                      basis$S, k)
    lp_with_grad <- find_loglikelihood_pen_with_grad(par, basis$X, data$y, data$c - 1, sp,
                                                     basis$S, k)

    lp_fun <- function(par) {
        find_loglikelihood_pen_noderiv(par, basis$X, data$y, data$c - 1, sp, basis$S, k)
    }

    library(numDeriv)
    grad_man_p <- grad(lp_fun, par)
    hess_man_p <- hessian(lp_fun, par)

    expect_equal(attr(lp_with_deriv, "gradient"), attr(lp_with_grad, "gradient"))
    expect_equal(attr(lp_with_deriv, "gradient"), grad_man_p)
    expect_equal(attr(lp_with_deriv, "hessian"), hess_man_p)

    library(microbenchmark)
    microbenchmark(
        l_pen_R(par, sp),
        find_loglikelihood_pen_noderiv(par, basis$X, data$y, data$c - 1, sp, basis$S, k),
        find_loglikelihood_pen_with_grad(par, basis$X, data$y, data$c - 1, sp, basis$S, k),
        find_loglikelihood_pen_with_hess(par, basis$X, data$y, data$c - 1, sp, basis$S, k)
    )
    #' timings similar to passing in beta: doing transform does not add too much cost
    #' could do ~ 1400 iterations per second with the gradient (when compiled outside of package)
    #' Now ~ 125 iterations per second with the gradient (when compiled in the package)
    #' Probably need to check compilation options to improve speed




})
