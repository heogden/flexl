test_that("Householder matrices work", {
    set.seed(1)
    
    u <- rnorm(10)
    H <- new_Householder(u)

    x <- rnorm(10)

    ## should just be able to just use as.matrix here
    ## but doesn't seem to be working
    H_mat <- as.matrix.Householder(H)
    H_mat_x <- as.numeric(H_mat %*% x)
    H_x <- find_Hx(H, x)

    expect_equal(H_mat_x, H_x)
    
})
