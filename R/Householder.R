#' not actually Householder matrices, but drop first column
#' to find H_star x such that H_star x = H(c(0, x)), where H is the
#' Householder matrix
new_Householder <- function(u) {
    stopifnot(is.numeric(u))
    gamma <- 2 / sum(u^2)
    structure(list(u = u, gamma = gamma), class = "Householder")
}

as.matrix.Householder <- function(x, ...) {
    u <- x$u
    gamma <- x$gamma

    H <- diag(nrow = length(u)) - gamma * outer(u, u)
    H[ , -1, drop = FALSE]
}

#' multiply H_star by x
find_Hstarx <- function(Hstar, x) {
    a <- sum(Hstar$u[-1] * x)
    c(0, x) - a * H_star$gamma * H_star$u
}
