drop_attributes <- function(x) {
    attributes(x) <- NULL
    x
}


log_det <- function(x) {
    drop_attributes(
        determinant(x, logarithm = TRUE)$modulus
    )
}

#' the log of the generalized determinant
log_det_gen <- function(x, rank) {
    evs <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    sum(log(evs[1:rank]))
}

is_neg_def <- function(A) {
    all(eigen(A, only.values = TRUE)$values < 0)
}
