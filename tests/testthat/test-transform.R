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

    dbeta_man <- array(dim = c(nbasis, length(alpha), k))
    for(i in 1:k) {
        dbeta_man[,,i] <- numDeriv::jacobian(function(x) {
            find_beta(x, nbasis, k)[,i]
        }, alpha)
    }
    
    dbeta_man_2 <- numDeriv::jacobian(find_beta, alpha, nbasis = nbasis, k = k)
    #' just binding together the rows of the three components of dbeta_man

    #' look at components
    alpha <- alpha[components == 1]
    Hstar <- find_Hstar(alpha)
    
    find_Hstar_x2 <- function(alpha, x2) {
        Hstar <- find_Hstar(alpha)
        
        a <- sum(Hstar$u[-1] * x2)
        c(0, x) - a * Hstar$gamma * Hstar$u
    }

    x <- rnorm(4)
    D_Hstar_x_man <- numDeriv::jacobian(find_Hstar_x2, alpha, x2 = x)
    D_Hstar_x <- attr(find_Hstar_x(Hstar, x, deriv = TRUE), "gradient")
    
    expect_equal(D_Hstar_x, D_Hstar_x_man)
    
}

