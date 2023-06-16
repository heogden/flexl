test_that("Transform gives orthogonal columns", {
    nbasis <- 5
    k <- 3

    components <- find_alpha_components(nbasis, k)

    
    set.seed(1)
    alpha <- rnorm(length(components))

    beta <- find_beta(alpha, nbasis, k)

    expect_equal(sum(beta[,2] * beta[,1]), 0)
    expect_equal(sum(beta[,3] * beta[,1]), 0)
    expect_equal(sum(beta[,3] * beta[,2]), 0)
})


test_that("can differentiate transform", {
    nbasis <- 5
    k <- 3

    components <- find_alpha_components(nbasis, k)

    set.seed(1)
    alpha <- rnorm(length(components))

    beta <- find_beta(alpha, nbasis, k)

    
    dbeta_man <- array(dim = c(nbasis, length(alpha), k))
    for(i in 1:k) {
        dbeta_man[,,i] <- numDeriv::jacobian(function(x) {
            find_beta(x, nbasis, k)$value[,i]
        }, alpha)
    }

    expect_equal(beta$gradient, dbeta_man)
}

