predict_flexl <- function(mod, newdata) {
    if(length(newdata$c) > 1)
        stop("can only predict for a single cluster at once")
    f0_x <- mod$f0(newdata$x)
    f_x <- mod$f(newdata$x)

    if(ncol(f_x) == 0)
        return(f0_x)
    
    row <- which(rownames(mod$u) == newdata$c)
    
    u <- mod$u[row,]
    f0_x + colSums(u * t(f_x))
}

