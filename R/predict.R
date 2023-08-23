predict_flexl <- function(mod, newdata) {
    if(length(newdata$c) > 1)
        stop("can only predict for a single cluster at once")

    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }

    newdata$x <- (newdata$x - mod$norm$m_x) / mod$norm$s_x
    
    f0_x <- mod$f0(newdata$x)

    if(mod$k == 0)
        return(f0_x)

    f_x <- mod$f(newdata$x)

    row <- which(rownames(mod$u) == newdata$c)
    
    u <- mod$u[row,]
    y_norm <- f0_x + colSums(u * t(f_x))
    mod$norm$m_y + mod$norm$s_y * y_norm
}

fitted_flexl <- function(mod, data) {
    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    
    u_ext <- mod$u[data$c, , drop = FALSE]
    y_norm <- unname(as.numeric(mod$f0_x) + rowSums(u_ext * mod$f_x))
    mod$norm$m_y + mod$norm$s_y * y_norm
}
