new_Householder <- function(u) {
    stopifnot(is.numeric(u))
    gamma <- 2 / sum(u^2)
    structure(list(u = u, gamma = gamma), class = "Householder")
}

as.matrix.Householder <- function(x, ...) {
    u <- x$u
    gamma <- x$gamma

    diag(nrow = length(u)) - gamma * outer(u, u)
}

#' multiply H by x
find_Hx <- function(H, x) {
    a <- sum(H$u * x)
    x - a * H$gamma * H$u
}

#' multiply H by c(0, x)
find_H0x <- function(H, x) {
    find_Hx(H, c(0, x))
}

                         
