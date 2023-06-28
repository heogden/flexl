predict_flexl <- function(mod, newdata) {
    if(length(newdata$c) > 1)
        stop("can only predict for a single cluster at once")
    f0_x <- mod$f0(newdata$x)

    if(mod$k == 0)
        return(f0_x)

    f_x <- mod$f(newdata$x)

    row <- which(rownames(mod$u) == newdata$c)
    
    u <- mod$u[row,]
    f0_x + colSums(u * t(f_x))
}

fitted_flexl <- function(mod, data) {
    u_ext <- mod$u[data$c, , drop = FALSE]
    unname(as.numeric(mod$f0_x) + rowSums(u_ext * mod$f_x))
}
