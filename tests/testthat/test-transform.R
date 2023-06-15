test_that("Transform gives orthogonal columns", {
    nbasis <- 5
    k <- 3

    components <- find_alpha_components(nbasis, k)

    
    set.seed(1)
    alpha <- rnorm(length(components))
    alpha_1 <- alpha[components == 1]
    alpha_2 <- alpha[components == 2]
    alpha_3 <- alpha[components == 3]


    beta_1 <- alpha_1
    
    H1 <- find_Hi(alpha_1)
    beta_2 <- find_Hx(H1, c(0, alpha_2))

    ## is beta_2 orthogonal to beta_1?
    expect_equal(sum(beta_1 * beta_2), 0)

    H2 <- find_Hi(alpha_2)
    beta_3 <- find_Hx(H1, c(0, find_Hx(H2, c(0, alpha_3))))

    expect_equal(sum(beta_1 * beta_3), 0)
    expect_equal(sum(beta_2 * beta_3), 0)
})


test_that("can differentiate transform", {
    
}

