find_orthogonal_complement_transform <- function(beta) {
    A <- qr(beta)
    k <- ncol(beta)
    qr.Q(A, complete = TRUE)[,-(1:k)]
}
