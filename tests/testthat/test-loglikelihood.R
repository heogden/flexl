test_that("derivatives of loglikelihood are correct", {
    data <- generate_test_data_1()$data
    
    nbasis <- 5
    sp <- 100
    basis <- find_orthogonal_spline_basis(nbasis, data$x)


    k <- 1

    row_list <- split(1:nrow(data), data$c)
    alpha_components <- find_alpha_components(nbasis, k)
    
    set.seed(1)
    par <- c(rnorm(nbasis + length(alpha_components), sd = 0.1), 1)

    par_split <- split_par(par, basis$nbasis)

    beta0 <- par_split$beta0
    beta <- find_beta(par_split$alpha, basis$nbasis, k)$value
    sigma <- exp(par_split$lsigma)

    f0 <- as.numeric(basis$X %*% beta0)
    fx <- basis$X %*% beta
    rows <- row_list[[1]]
    y <- data$y
    
    l1 <- find_loglikelihood_cluster(rows, f0, fx, y, sigma)
    l2 <- find_loglikelihood_cluster_struc(rows, f0, fx, y, sigma)
    l3 <- find_loglikelihood_cluster_from_Cpp(rows, f0, fx, y, sigma)
    expect_equal(l2, l1)
    expect_equal(l3, l1)

    library(microbenchmark)
    microbenchmark(
        find_loglikelihood_cluster(rows, f0, fx, y, sigma),
        find_loglikelihood_cluster_struc(rows, f0, fx, y, sigma),
        find_loglikelihood_cluster_from_Cpp(rows, f0, fx, y, sigma)
    )
    #' new method in R is a little bit slower for k > 1, faster for k = 1
    #' new method in C++ is the fastest of all

    
    X <- basis$X
    X_c <- X[rows, ]
    y_c <- y[rows]

    l4 <- find_loglikelihood_cluster_Stan(par, X_c, y_c)

    l_cluster_basic <- function(par) {
        beta_0 <- par[1:5]
        beta_1 <- par[6:10]
        sigma <- par[11]
        f0_c <- X_c %*% beta_0
        f1_c <- X_c %*% beta_1
        Sigma <- diag(sigma^2, nrow = 10) + outer(as.numeric(f1_c), as.numeric(f1_c))
        mvtnorm::dmvnorm(y_c, f0_c, Sigma, log = TRUE)
    }

    library(numDeriv)
    grad_man <- grad(l_cluster_basic, par)
    hess_man <- hessian(l_cluster_basic, par)

    expect_equal(attr(l4, "gradient"), grad_man)
    expect_equal(attr(l4, "hessian"), hess_man)

    
    #' what if cluster size was much larger?
    #' e.g. for whole data?

    rows <- Reduce(c, row_list)

    microbenchmark(
        find_loglikelihood_cluster(rows, f0, fx, y, sigma),
        find_loglikelihood_cluster_struc(rows, f0, fx, y, sigma),
        find_loglikelihood_cluster_from_Cpp(rows, f0, fx, y, sigma)
    )
    #' if we have a larger number of obs/cluster, new method in R wins even for larger k
    #' for some reason, new method in C++ is slightly slower than in R for large number
    #' of obs/cluster: because we are passing the vectors/matrices about?
    
    pen_deviance_grad_man <- numDeriv::grad(find_pen_deviance, par, sp = sp, y = data$y,
                                            row_list = row_list, basis = basis, k = k)
    pen_deviance_grad_auto <- attr(pen_deviance, "gradient")
    expect_equal(pen_deviance_grad_auto, pen_deviance_grad_man)

    pen_deviance_hess_man <- numDeriv::hessian(find_pen_deviance, par, sp = sp, y = data$y,
                                               row_list = row_list, basis = basis, k = k)
    pen_deviance_hess_auto <- attr(pen_deviance, "hessian")
    expect_equal(pen_deviance_hess_auto, pen_deviance_hess_man)


})

test_that("can use Stan with autodiff", {
    set.seed(1)
    n <- 100
    x <- rnorm(n, mean = 0, sd = 1)

    ldnorm_with_derivs(x, c(1, 1))
})

