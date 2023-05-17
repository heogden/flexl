test_that("can differentiate transform", {
    nbasis <- 5

    set.seed(1)
    beta <- matrix(rnorm(nbasis), ncol = 1)

    T <- find_orthogonal_complement_transform(beta)

    A <- qr(beta)
    k <- ncol(beta)
    Q <- qr.Q(A, complete = TRUE)
    T <- Q[,-(1:k)]
    #' Can I find this QR decomposition manually,
    #' by using Householder reflections?

    #' From Wood2017, B.6, need to find Householder matrix H1
    #' (then Q = H1 in this case, since m = 1)
    #' To find H1, set u = beta - (||beta||, 0, .., 0)^T

    find_H1 <- function(beta) {
        beta_norm <- sqrt(sum(beta^2))
        n <- nrow(beta)
        e1 <- c(1, rep(0, n-1))
        u <- beta - beta_norm * e1
        gamma <- 2 / sum(u^2)
        diag(nrow = n, ncol = n) - gamma * tcrossprod(u)
    }
    

    H1 <- find_H1(beta)
    expect_equal(H1, Q)


    find_dH1 <- function(beta) {
        beta_norm <- sqrt(sum(beta^2))
        n <- nrow(beta)
        e1 <- c(1, rep(0, n-1))
        u <- beta - beta_norm * e1
        gamma <- 2 / sum(u^2)
        du <- diag(nrow = n, ncol = n)
        beta_norm <- sqrt(sum(beta^2))
        du[1,] <- du[1,] - beta/beta_norm
        dgamma <- as.numeric(-gamma^2 * crossprod(du, u))

        uut <- tcrossprod(u)
        u <- as.numeric(u)
        -outer(uut, dgamma) - gamma * (outer(u, du) + aperm(outer(u, du), c(2, 1, 3)))
    }

    dH1_man <- numDeriv::jacobian(find_H1, beta)
    dH1 <- find_dH1(beta)
    
    expect_equal(as.numeric(dH1), as.numeric(dH1_man))

   
    dT_man <- numDeriv::jacobian(find_orthogonal_complement_transform, beta)
    dT <- dH1[, -1, ]
    expect_equal(as.numeric(dT), as.numeric(dT_man))
 
})
