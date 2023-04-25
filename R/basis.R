find_orthogonal_spline_basis <- function(nbasis, x) {
    interior_knots <- seq(min(x) - 0.01 * sd(x), max(x) + 0.01 * sd(x), length = nbasis - 2)
    knots <- orthogonalsplinebasis::expand.knots(interior_knots)
    basis <- orthogonalsplinebasis::OrthogonalSplineBasis(knots)
    list(X = orthogonalsplinebasis::evaluate(basis, x = x),
         S = orthogonalsplinebasis::OuterProdSecondDerivative(basis),
         basis = basis)
}

#' if x outside knots, extrapolate constant at nearest knot point, rather than NA
#' (which is the default in orthogonalsplinebasis:::evaluate.SplineBasis
evaluate_with_extrap <- function(object, x) {
    knots <- object@knots
    order <- object@order
    lower <- knots[order]
    upper <- knots[length(knots)-order+1]
    x[x<lower] <- lower
    x[x>upper] <- upper
    orthogonalsplinebasis::evaluate(object, x)
}


find_spline_fun <- function(beta, basis) {
    spline_fun <- function(x) {
        B <-  evaluate_with_extrap(basis$basis, x = x)
        result <- tcrossprod(B, t(beta))
        if(!is.matrix(beta))
            result <- as.numeric(result)
        result
    }
    spline_fun
}

