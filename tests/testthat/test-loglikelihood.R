test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    k <- 2
    sp <- 10
    basis <- find_orthogonal_spline_basis(nbasis, data$x)

    alpha_components <- find_alpha_components(nbasis, k)

    set.seed(1)
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), 1)

    lp_grad <- loglikelihood_pen_grad(par, basis$X, data$y, data$c - 1, sp, basis$S, k)
    lp_hess <- loglikelihood_pen_hess(par, basis$X, data$y, data$c - 1, sp, basis$S, k)

    lp_fun <- function(par) {
        loglikelihood_pen(par, basis$X, data$y, data$c - 1, sp, basis$S, k)
    }

    library(numDeriv)
    lp_grad_man <- grad(lp_fun, par)
    lp_hess_man <- hessian(lp_fun, par)

    expect_equal(lp_grad, lp_grad_man)
    expect_equal(lp_hess, lp_hess_man)

    library(microbenchmark)
    microbenchmark(
        loglikelihood_pen(par, basis$X, data$y, data$c - 1, sp, basis$S, k),
        loglikelihood_pen_grad(par, basis$X, data$y, data$c - 1, sp, basis$S, k),
        loglikelihood_pen_hess(par, basis$X, data$y, data$c - 1, sp, basis$S, k)
    )
    #' timings similar to passing in beta: doing transform does not add too much cost
    #' can do ~ 900 iterations per second with the gradient



})
