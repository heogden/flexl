find_orthogonal_complement_transform <- function(beta) {
    A <- qr(beta)
    k <- ncol(beta)
    qr.Q(A, complete = TRUE)[,-(1:k)]
}

#' Code below modified from
#' Differentiating the QR Decomposition
#' Jan de Leeuw
#' May 2023
#' DOI:10.13140/RG.2.2.19201.12640
lt <- function(x) {
    n <- nrow(x)
    x[outer(1:n, 1:n, "<")] <- 0
    return(x)
}

d_qr <- function(x, y) {
    n <- nrow(x)
    m <- ncol(x)
    z <- x + y
    qrx <- qr(x)
    qx <- qr.Q(qrx)
    rx <- qr.R(qrx)
    qrz <- qr(z)
    qz <- qr.Q(qrz)
    rz <- qr.R(qrz)
    qp <- qr.Q(qr(cbind(qx, diag(n))))[, -(1:m)]
    ri <- solve(rx)
    v <- crossprod(qx, y %*% ri)
    a <- lt(v) - t(lt(v))
    b <- crossprod(qp, y %*% ri)
    dq <- qx %*% a + qp %*% b
    dr <- (v - a) %*% rx
    return(list(
        qx = qx,
        rx = rx,
        qz = qz,
        rz = rz,
        dq = dq,
        dr = dr
    ))
}

p_q <- function(x) {
    n <- nrow(x)
    m <- ncol(x)
    qrx <- qr(x)
    qx <- qr.Q(qrx)
    rx <- qr.R(qrx)
    ri <- solve(rx)
    qp <- qr.Q(qr(cbind(qx, diag(n))))[, -(1:m), drop = FALSE]
    g <- array(NA, dim = c(n, m, n, m))
    for (i in 1:n) {
        for (j in 1:m) {
            v <- outer(qx[i, ], ri[j, ])
            a <- lt(v) - t(lt(v))
            b <- outer(qp[i, ], ri[j, ])
            dq <- qx %*% a + qp %*% b
            ##dr <- (v - a) %*% rx
            ##g[, k] <- c(as.vector(dq), as.vector(dr))
            g[, , i, j] <- dq
        }
    }
    return(g)
}


#' From https://discourse.julialang.org/t/using-forwarddiff-on-qr-factorization/42625/12
## d_q <- function(A) {
##     m <- nrow(A)
##     n <- ncol(A)
##     QR = qr(A)
##     Q <- qr.Q(QR)
##     QdA = Q'*dA
##     QdAR = QdA/R
##     X = tril(QdAR,-1)
##     X = X - X'
##     dQ = Q*X + dA/R - Q*QdAR
##     dR = QdA - X*R
##     return dQ,dR
## end
