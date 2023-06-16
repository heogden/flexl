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
    
    

    find_Hstar_x <- function(alpha, x2) {
        alpha_norm <- sqrt(sum(alpha^2))
    
        t <- alpha_norm * (alpha_norm - alpha[1])
        dt <- 2 * (1 - alpha[1]) / alpha_norm * alpha
        dt[1] <- dt[1] - alpha_norm
        
        gamma <- 1 / t
        dgamma <- - dt / t^2
        
        u <- alpha
        u[1] <- u[1] - alpha_norm
        a <- sum(u[-1] * x2)

        c(0, x2) - a * gamma * u
    }

    x <- rnorm(4)
    D_Hstar_x_man <- numDeriv::jacobian(find_Hstar_x, alpha, x2 = x)

    
    find_D_Hstar_x <- function(alpha, x) {
         alpha_norm <- sqrt(sum(alpha^2))

         u <- alpha
         u[1] <- u[1] - alpha_norm

         t <- sum(u^2)
         dt <- 4 * alpha - 2 * alpha[1] * alpha / alpha_norm
         dt[1] <- dt[1] - 2 * alpha_norm
         
         gamma <- 2 / t
         dgamma <- - 2 * dt / t^2
         
         a <- sum(u[-1] * x)

         n <- length(alpha)
         du <- diag(nrow = n)
         du[1,] <- du[1,] - alpha / alpha_norm

         M <- outer(u, dgamma) + gamma * du

         - gamma * outer(u, c(0, x)) - a * M
    }

    D_Hstar_x <- find_D_Hstar_x(alpha, x)
    expect_equal(D_Hstar_x, D_Hstar_x_man)


    
}

