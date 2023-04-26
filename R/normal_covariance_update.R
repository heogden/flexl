#' Find log phi(z; 0, Sigma_k),
#' where z is fixed and stored
#' and
#' Sigma_k = Sigma_km1 + a a^T
#' for each new value of a.
#' Here info_km1 contains values needed corresponding to Sigma_km1 and z
ldmvnorm <- function(a, info_km1) {
    b <- info_km1$Sigma_inv %*% a
    c <- 1 + sum(a * b)
    ldet_Sigma <- log(c) + info_km1$ldet_Sigma
    d <- sum(info_km1$Sigma_inv_z * a)

    n <- length(info_km1$z)

    quad_form <- info_km1$tz_Sigma_inv_z - d^2/c
    -n/2 * log(2*pi) - 1/2 * ldet_Sigma - 1/2 * quad_form
}

#' Derivative of log(|Sigma_k|) with respect to a
#' (computed above as const + log c)
log_det_Sigma_grad <- function(c, dc) {
    dc / c
}

#' Derivative of z^T (Sigma_k)^{-1} z with respect to a
#' (computed above as const - d^2 / c)
quad_form_grad <- function(c, dc, d, dd) {    
    - 2 * d * dd / c + d^2 * dc / c^2
}

#' could do this at the same time as computing values
#' and avoid recomputing quantities
ldmvnorm_grad <- function(a, info_km1) {
    b <- as.numeric(info_km1$Sigma_inv %*% a)
    c <- 1 + sum(a * b)
    dc <- 2 * b
    d <- sum(info_km1$Sigma_inv_z * a)
    dd <- as.numeric(info_km1$Sigma_inv_z)

    - 0.5 * (log_det_Sigma_grad(c, dc) + quad_form_grad(c, dc, d, dd))
}

log_det_Sigma_hess <- function(c, dc, d2c) {
    d2c / c - tcrossprod(dc / c)
}

quad_form_hess <- function(c, dc, d2c, d, dd) {
    t1 <- -2 * tcrossprod(dd) / c
    t2 <- 2 * d * tcrossprod(dd, dc) / c^2
    t3 <- 2 * d * tcrossprod(dc, dd) / c^2
    t4 <- d^2 * d2c / c^2
    t5 <- -2 * d^2 * tcrossprod(dc) / c^3

    t1 + t2 + t3 + t4 + t5
}


ldmvnorm_hess <- function(a, info_km1) {
    b <- as.numeric(info_km1$Sigma_inv %*% a)
    
    c <- 1 + sum(a * b)
    dc <- 2 * b
    d2c <- 2 * info_km1$Sigma_inv
    
    d <- sum(info_km1$Sigma_inv_z * a)
    dd <- as.numeric(info_km1$Sigma_inv_z)

    -0.5 * (log_det_Sigma_hess(c, dc, d2c) + quad_form_hess(c, dc, d2c, d, dd))
}
