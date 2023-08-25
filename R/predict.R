predict_y_given_mod <- function(mod, newdata) {
    newdata_norm <- newdata
    newdata_norm$x <- (newdata$x - mod$norm$m_x) / mod$norm$s_x

    x_norm_all <- sort(unique(newdata_norm$x))
    newdata_norm_all <- tidyr::crossing(c = unique(mod$data$c), x = x_norm_all)
    y_norm_fun <- find_spline_fun(mod$par_cluster, mod$basis)
    ynorm_all <- y_norm_fun(x_norm_all)
    newdata_norm_all$y_hat <- as.numeric(ynorm_all)

    y_norm <- merge(newdata_norm, newdata_norm_all, by = c("c", "x"), sort = FALSE)$y_hat
    
    mod$norm$m_y + mod$norm$s_y * y_norm
}

predict_y_given_sample <- function(sample, mod, newdata) {
    mod_s <- mod
    mod_s$par_cluster <- sample
    
    predict_y_given_mod(mod_s, newdata)
}



predict_flexl <- function(mod, newdata, interval = "none",
                          level = 0.95, samples = NULL,
                          n_samples = 1000) {
    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    y_hat <- predict_y_given_mod(mod, newdata)
    if(interval == "none") {
        return(y_hat)
    } else {
        if(is.null(samples)) {
            samples <- find_samples(mod, n_samples)
        }
        y_hat_samples <- sapply(samples, predict_y_given_sample,
                                mod = mod, newdata = newdata)
        alpha <- (1 - level) / 2
        y_hat_lower <- apply(y_hat_samples, 1, quantile, probs = alpha)
        y_hat_upper <- apply(y_hat_samples, 1, quantile, probs = 1 - alpha)
        return(data.frame(estimate = y_hat, lower = y_hat_lower, upper = y_hat_upper))
    }
    
   
}


fitted_flexl <- function(mod) {
    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    
    u_ext <- mod$u_hat[mod$data$c, , drop = FALSE]
    y_norm <- unname(as.numeric(mod$f0_x) + rowSums(u_ext * mod$f_x))
    mod$norm$m_y + mod$norm$s_y * y_norm
}
