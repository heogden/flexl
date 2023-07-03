test_that("compute log density from eigenvalues correctly", {
    n <- 20
    x <- seq(0, 1, length.out = n)
    f1 <- rep(1, n)
    f2 <- x - 0.5 * mean(f1) #' to make orthogonal
    sigma2 <- 0.5
    Sigma <- outer(f1, f1) + outer(f2, f2) + diag(sigma2, n, n)


    Sigma_inv <- solve(Sigma)
    tau <- 1 / sigma2
    Sigma_inv_rem <- Sigma_inv - diag(tau, n, n)
    d <- sqrt(-diag(Sigma_inv_rem))
    Sigma_inv_rem_pred <- -outer(d, d)
    expect_equal(-outer(d, d), Sigma_inv_rem)

    which.max(abs(-outer(d, d) - Sigma_inv_rem))

    
    
    evalues <- eigen(Sigma)$values

    #' manually finding the eigenvalues
    #' (only works if fk are orthogonal)
    lambda1 <- sum(f1^2) + sigma2
    lambda2 <- sum(f2^2) + sigma2

    evalues_auto <- c(lambda1, lambda2, rep(sigma2, n-2))
    expect_equal(evalues_auto, evalues)
    
    evectors <- eigen(Sigma)$vectors

    
    #' manually finding the eigenvectors
    evectors[,1]
    evectors[,1] / f1 #' a constant, chosen to make f1 normal
    sum(evectors[,1]^2)
    1 / sqrt(sum(f1^2)) #' gets constant (up to sign)

    evectors[,2]
    evectors[,2] / f2 #' a constant
    sum(evectors[,2]^2)
    1 / sqrt(sum(f2^2)) #' gets constant


    #' possibly rest of eigenvectors fairly arbitrary?
    #' but they do need to be orthogonal to first two eigenvectors
    #' Could be chosen multiple ways?
    evectors[,3]
    evectors[,4]

    library(mvtnorm)
    set.seed(1)
    z <- as.numeric(rmvnorm(1, sigma = Sigma))
    ld <- dmvnorm(z, sigma = Sigma, log = TRUE)
    ld_e <- ldmvnorm_eigen(z, evalues)
    expect_equal(ld_e, ld)

    #' look at parts separately
    evalues_inv <- 1 / evalues
    
    Sigma_inv <- solve(Sigma)
    const <- - n / 2 * log(2 * pi)
    
    log_det_Sigma_inv <- determinant(Sigma_inv)$modulus
    log_det_Sigma_inv_e <- sum(log(evalues_inv))
    #' these are the same
    
    zT_Sigma_inv_z <- as.numeric(t(z) %*% Sigma_inv %*% z)
    zT_Sigma_inv_z_e <- sum(z^2 * evalues_inv)
    #' this is the problem part

    
    const + log_det_Sigma_inv / 2 - zT_Sigma_inv_z / 2
    
})
