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

    find_beta_from_alpha <- function(alpha, nbasis, k) {
        alpha_split <- split_alpha(alpha, nbasis, k)
        find_beta(alpha_split)
    }
    set.seed(1)
    alpha <- rnorm(3 * nbasis - 1)
    k <- 2
    beta <- find_beta_from_alpha(alpha, nbasis, k)

    #' this is what we want to end up with
    #' (ignore first nbasis elements, giving alpha_0, they just control beta_0)
    dbeta_man <- numDeriv::jacobian(find_beta_from_alpha, alpha, nbasis = nbasis, k = k)[,-(1:nbasis)]
    #' beta_1 influenced by alpha_1 only (in fact beta_1 = alpha_1). So deriv is identity
    #' beta_2 influenced by alpha_1 and alpha_2
    #' Focus on this part derivatives of beta_2 write alpha = (alpha_1, alpha_2)
    dbeta2_man <- dbeta_man[-(1:nbasis), ]

    T <- find_orthogonal_complement_transform(beta[, 1, drop = FALSE])
    alpha_12 <- alpha[-(1:nbasis)]
    alpha_1 <- alpha_12[1:nbasis]
    alpha_2 <- alpha_12[-(1:nbasis)]
    beta_2 <- T %*% alpha_2

    #' derivative of beta_2 wrt alpha_2 is just T
    #' so what is left is derivatives of beta_2 wrt alpha_1
    #' T is a function of beta_1 = alpha_1
    dH1 <- find_dH1(beta[, 1, drop = FALSE])
    dT <- dH1[, -1, ]
    
    dbeta2_alpha1 <- matrix(nrow = nbasis, ncol = nbasis)
    for(i in 1:nbasis)
        for(j in 1:nbasis)
            dbeta2_alpha1[i,j] <- sum(dT[i,,j] * alpha_2)

    dbeta2_alpha1_man <- dbeta2_man[,1:nbasis]
    expect_equal(dbeta2_alpha1, dbeta2_alpha1_man)
})
