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
quad_form_grad <- function(c, d, dc, dd) {    
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

    - 0.5 *(log_det_Sigma_grad(c, dc) + quad_form_grad(c, d, dc, dd))
}
