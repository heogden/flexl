predict_y_given_mod <- function(mod, newdata) {
    newdata$x <- (newdata$x - mod$norm$m_x) / mod$norm$s_x
    
    f0_x <- mod$f0(newdata$x)

    if(mod$k == 0)
        return(f0_x)

    f_x <- mod$f(newdata$x)

    row <- which(rownames(mod$u_hat) == newdata$c)
    
    u <- mod$u_hat[row,]
    y_norm <- f0_x + colSums(u * t(f_x))
    mod$norm$m_y + mod$norm$s_y * y_norm
}


predict_y_given_sample <- function(sample, mod, newdata) {
    mod_s <- sample
    mod_s$k <- mod$k
    mod_s$norm <- mod$norm
    
    predict_y_given_mod(mod_s, newdata)
}




predict_flexl <- function(mod, newdata, interval = "none",
                          level = 0.95, samples = NULL,
                          n_samples = 1000) {
    if(length(newdata$c) > 1)
        stop("can only predict for a single cluster at once")

    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    if(interval == "none") {
        predict_y_given_mod(mod, newdata)
    } else {
        if(is.null(samples)) {
            samples <- find_samples(mod, n_samples)
        }
        y_hat_samples <- sapply(samples, predict_y_given_sample,
                                mod = mod, newdata = newdata)
        
    }
    
    

   
}

fitted_flexl <- function(mod, data) {
    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    
    u_ext <- mod$u[data$c, , drop = FALSE]
    y_norm <- unname(as.numeric(mod$f0_x) + rowSums(u_ext * mod$f_x))
    mod$norm$m_y + mod$norm$s_y * y_norm
}
