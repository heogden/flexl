test_that("Log ML not always increasing with k", {
    
    find_log_ml <- function(sp, k, data, nbasis) {
        basis <- find_orthogonal_spline_basis(nbasis, data$x)
        S <- basis$S
        
        mod <- fit_given_sp_init(data, sp, k, basis)
        mod$log_ml
    }

 
    nbasis <- 10
    
    #' generate data with truth k = 0
    data0 <- generate_test_data_0()$data

    sp_poss <- seq(0.01, 1, length.out = 100)
    log_ml_0_poss <- sapply(sp_poss, find_log_ml, k = 0, data = data0, nbasis = nbasis)
    log_ml_1_poss <- sapply(sp_poss, find_log_ml, k = 1, data = data0, nbasis = nbasis)
    log_ml_2_poss <- sapply(sp_poss, find_log_ml, k = 2, data = data0, nbasis = nbasis)

    plot(range(sp_poss), range(c(log_ml_0_poss, log_ml_1_poss, log_ml_2_poss)), type = "n")
    lines(sp_poss, log_ml_0_poss)
    lines(sp_poss, log_ml_1_poss, lty = 2)
    lines(sp_poss, log_ml_2_poss, lty = 3)
        

    expect_gt(max(log_ml_0_poss), max(log_ml_1_poss))
    expect_gt(max(log_ml_1_poss), max(log_ml_2_poss))

    
    #' generate data with truth k = 2
    data2 <- generate_test_data_2()$data
    
    sp_poss <- seq(0.01, 10, length.out = 100)
    log_ml_0_poss <- sapply(sp_poss, find_log_ml, k = 0, data = data2, nbasis = nbasis)
    log_ml_1_poss <- sapply(sp_poss, find_log_ml, k = 1, data = data2, nbasis = nbasis)
    log_ml_2_poss <- sapply(sp_poss, find_log_ml, k = 2, data = data2, nbasis = nbasis)
    log_ml_3_poss <- sapply(sp_poss, find_log_ml, k = 3, data = data2, nbasis = nbasis)

    plot(range(sp_poss),
         range(c(log_ml_0_poss, log_ml_1_poss, log_ml_2_poss, log_ml_3_poss)), type = "n")
    lines(sp_poss, log_ml_0_poss)
    lines(sp_poss, log_ml_1_poss, lty = 2)
    lines(sp_poss, log_ml_2_poss, lty = 3)
    lines(sp_poss, log_ml_3_poss, lty = 4)


    expect_gt(max(log_ml_1_poss), max(log_ml_0_poss))
    expect_gt(max(log_ml_2_poss), max(log_ml_1_poss))
    expect_gt(max(log_ml_2_poss), max(log_ml_3_poss))


})
