test_that("find norm of function from basis coefficients", {
   nbasis <- 5

   n <- 100
   x <- runif(n)

   basis <- find_orthogonal_spline_basis(nbasis, x)

   beta <- matrix(rnorm(nbasis), ncol = 1)
   f <- find_spline_fun(beta, basis)

   f_norm2 <- sum(beta^2)
   f_norm2_man <- integrate(function(x){f(x)^2}, min(basis$basis@knots), max(basis$basis@knots))$value

   expect_equal(f_norm2, f_norm2_man)
   
})

test_that("find deriv of function written in spline basis", {
    nbasis <- 5

    n <- 100
    x <- runif(n)
    
    basis <- find_orthogonal_spline_basis(nbasis, x)
    
    beta <- matrix(rnorm(nbasis), ncol = 1)
    f <- find_spline_fun(beta, basis)
    
    x_test <- c(-0.5, 0.2, 0.3)

    grad_f_man <- numDeriv::grad(f, x_test)
    grad_f <- f(x_test, deriv = TRUE)

    expect_equal(as.numeric(grad_f), grad_f_man)
}

