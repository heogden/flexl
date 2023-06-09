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
    k <- 3
    n_alpha_each <- c(nbasis, nbasis - 0:(k-1))
    alpha <- rnorm(sum(n_alpha_each))
    beta <- find_beta_from_alpha(alpha, nbasis, k)


    A <- qr(beta)
    k <- ncol(beta)
    Q <- qr.Q(A, complete = TRUE)
    R <- qr.R(A, complete = TRUE)
    T <- Q[,-(1:k)]
    
    #' this is what we want to end up with
    dbeta_man <- numDeriv::jacobian(find_beta_from_alpha, alpha, nbasis = nbasis, k = k)
    #' (can ignore first nbasis elements, giving alpha_0, they just control beta_0)
    #' beta_1 influenced by alpha_1 only (in fact beta_1 = alpha_1). So deriv is identity
    #' beta_2 influenced by alpha_1 and alpha_2
    #' Focus on this part derivatives of beta_2 write alpha = (alpha_1, alpha_2)
    beta_id <- rep(1:k, each = nbasis)
    alpha_id <- rep(0:k, n_alpha_each)
    dbeta2_man <- dbeta_man[beta_id == 2, alpha_id <= 2]

    alpha_1 <- alpha[alpha_id == 1]
    beta_1 <- alpha_1
    T <- find_orthogonal_complement_transform(matrix(beta_1, ncol = 1))

    alpha_2 <- alpha[alpha_id == 2]
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

    dbeta2_alpha1_man <- dbeta_man[beta_id == 2, alpha_id == 1]
    expect_equal(dbeta2_alpha1, dbeta2_alpha1_man)

    #' Now try to do this directly, based on
    #' B.5 of Wood2017:
    #' f1(beta_1, x) = H1 x,
    #' H1 is Householder transform of beta_1.
    f1 <- function(beta_1, x) {
        beta_1_norm <- sqrt(sum(beta_1^2))
        u <- beta_1
        u[1] <- beta_1[1] - beta_1_norm
        gamma <- 2 / sum(u^2)
        delta <- sum(u * x)
        x - delta * gamma * u
    }
    expect_equal(f1(alpha_1, c(0, alpha_2)), as.numeric(beta_2))

    #' To calculate derivatives of f1 with respect to beta_1
    #' could store u, gamma, delta and derivs,
    #' then can calculate f1 and deriv wrt beta_1 for different x quickly.
    df1 <- function(beta_1, x) {
        beta_1_norm <- sqrt(sum(beta_1^2))
        u <- beta_1
        u[1] <- beta_1[1] - beta_1_norm
        gamma <- 2 / sum(u^2)
        delta <- sum(u * x)

        n <- length(beta_1)
        du <- diag(nrow = n, ncol = n)
        du[1,] <- du[1,] - beta_1/beta_1_norm
        dgamma <- as.numeric(-gamma^2 * crossprod(du, u))
        ddelta <- crossprod(du, x)
        
        - gamma * tcrossprod(u, ddelta) - delta * tcrossprod(u, dgamma) - delta * gamma * du
    }

    beta_1 <- alpha_1
    expect_equal(df1(beta_1, c(0, alpha_2)), dbeta2_alpha1)

    dbeta3_alpha_2_man <- dbeta_man[beta_id == 3, alpha_id == 2]
    dbeta3_alpha_1_man <- dbeta_man[beta_id == 3, alpha_id == 1]
               
    alpha_3 <- alpha[alpha_id == 3]

    f2 <- function(beta_1, beta_2, x) {
        x2 <- f1(beta_1, x)
        c(x2[1], f1(beta_2[-1], x2[-1]))
    }
    
    
    f2(beta[,1], beta[,2], c(0, 0, alpha_3))
    beta[,3]
    #' not yet working

    H_12 <- find_H1(beta[-1, 2, drop = FALSE])
    H_2 <- rbind(0, cbind(0, H_12))
    H_2[1,1] <- 1

    H_1 <- find_H1(beta[, 1, drop = FALSE])

    as.numeric(H_2 %*% (H_1 %*% c(0, 0, alpha_3)))
    #' same as f_2(...), but not same as beta[,3]

    n <- length(beta_1)
    H_13 <- find_H1(beta[-c(1,2), 3, drop = FALSE])
    H_3 <- diag(nrow = n, ncol = n)
    H_3[-c(1,2), -c(1, 2)] <- H_13

    H_3 %*% (H_2 %*% (H_1 %*% beta))
    #' does not have correct upper triangular form
    #' something is going wrong here

    step1 <- H_1 %*% beta
    #' meets requirement that have zeroes on first col, except first element
    step2 <- H_2 %*% step1
    #' does **not** meet the requirement that second element of second column unchanged
    #' and everything beyond the second element is zeroed.

    #' create H2 based on step1 instead of beta?
    H_12 <- find_H1(step1[-1, 2, drop = FALSE])
    H_2 <- rbind(0, cbind(0, H_12))
    H_2[1,1] <- 1

    step2 <- H_2 %*% step1
    #' now sets everything beyond second element of second column to zero
    #' amd the first row is unchanged

    expect_equal(step1[,2], c(0, alpha_2))

    
    H_13 <- find_H1(step2[-c(1, 2), 3, drop = FALSE])
    H_3 <- rbind(0, 0, cbind(0, 0, H_13))
    H_3[1, 1] <- 1
    H_3[2, 2] <- 1

    step3 <- H_3 %*% step2

    Q_my <- H_1 %*% H_2 %*% H_3
    R_my <- step3
    R_my_nonzero <- R_my[1:ncol(R_my), ]

    expect_equal(Q_my %*% R_my, beta)
    expect_equal(crossprod(Q_my), diag(nrow = nrow(beta))) # Q is orthogonal
    expect_true(all(abs(R_my_nonzero[row(R_my_nonzero) < col(R_my_nonzero)]) < 1e-8))
                                        # Upper part of R is upper triangular

    
    T_my <- Q_my[,-(1:k)]
    ## T is different from T_my: we have a different transformation,
    ## both come from QR decompositions, but there are multiple possible choices

    k <- 3
    n_alpha_each <- c(nbasis, nbasis - 0:(k-1))
    alpha <- rnorm(sum(n_alpha_each))
    beta3 <- find_beta_from_alpha(alpha, nbasis, k)
    beta2 <- beta3[,1:2]

    ## Idea, we form transformation by repeated finding of Householder reflections
    ## Can we take H_1 and H_2 from beta2 in use to find tranform for beta3?
    ## Yes: because H_1 depends only on first column of beta, same for beta2 and beta3
}

