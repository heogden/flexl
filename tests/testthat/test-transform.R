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

    
    T_ij <- function(beta) {
        find_orthogonal_complement_transform(beta)[i, j]
    }

    
    dT_man <- array(NA, dim = c(dim(T), nbasis))
    
    for(i in 1:nrow(T))
        for(j in 1:ncol(T))
            dT_man[i, j, ] <- numDeriv::grad(T_ij, beta)

    #' can we recover this automatically?
    numDeriv::jacobian(find_orthogonal_complement_transform, beta)
    
    d_qr_out <- p_q(beta)

    d_qr_out[,,,1]
})
