find_orthogonal_spline_basis <- function(nbasis, x) {
    interior_knots <- seq(min(x) - 0.01 * sd(x), max(x) + 0.01 * sd(x), length = nbasis - 2)
    knots <- orthogonalsplinebasis::expand.knots(interior_knots)
    basis <- orthogonalsplinebasis::OrthogonalSplineBasis(knots)
    list(X = orthogonalsplinebasis::evaluate(basis, x = x),
         S = orthogonalsplinebasis::OuterProdSecondDerivative(basis$osb),
         basis = basis)
}
