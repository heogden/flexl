predict_y_given_mod <- function(mod, newdata, deriv = FALSE) {
    newdata_norm <- newdata
    newdata_norm$x <- (newdata$x - mod$norm$m_x) / mod$norm$s_x

    x_norm_all <- sort(unique(newdata_norm$x))
    clusters <- sort(unique(mod$data$c))
    newdata_norm_all <- tidyr::expand_grid(c = clusters, x = x_norm_all)
    y_norm_fun <- find_spline_fun(mod$par_cluster, mod$basis)
    ynorm_all <- y_norm_fun(x_norm_all, deriv)
    newdata_norm_all$y_hat <- as.numeric(ynorm_all)

    newdata_norm$id <- 1:length(newdata_norm$x)

    newdata_merged <- merge(newdata_norm, newdata_norm_all, by = c("c", "x"), sort = FALSE)
    y_norm <- newdata_merged$y_hat[order(newdata_merged$id)]

    if(deriv) {
        a <- 0
        b <- mod$norm$s_y / mod$norm$s_x
    }
    
    else {
        a <- mod$norm$m_y
        b <- mod$norm$s_y
    }
    
    a + b * y_norm
}

#' @export
predict_y_given_sample <- function(sample, mod, newdata, deriv = FALSE) {
    mod_s <- mod
    mod_s$par_cluster <- sample
    
    predict_y_given_mod(mod_s, newdata, deriv)
}


#' @export
predict_flexl <- function(mod, newdata, deriv = FALSE, interval = FALSE,
                          level = 0.95, samples = NULL,
                          n_samples = 1000) {
    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    y_hat <- predict_y_given_mod(mod, newdata, deriv)
    if(!interval) {
        return(y_hat)
    } else {
        if(is.null(samples)) {
            samples <- find_samples(mod, n_samples)
        }
        y_hat_samples <- sapply(samples, predict_y_given_sample,
                                mod = mod, newdata = newdata, deriv = deriv)
        alpha <- (1 - level) / 2
        y_hat_lower <- apply(y_hat_samples, 1, quantile, probs = alpha)
        y_hat_upper <- apply(y_hat_samples, 1, quantile, probs = 1 - alpha)
        return(data.frame(estimate = y_hat, lower = y_hat_lower, upper = y_hat_upper))
    }
    
   
}

#' @export
fitted_flexl <- function(mod) {
    if(is.null(mod$norm)) {
        mod$norm <- list(m_y = 0, s_y = 1, m_x = 0, s_x = 1)
    }
    
    u_ext <- mod$u_hat[mod$data$c, , drop = FALSE]
    y_norm <- unname(as.numeric(mod$f0_x) + rowSums(u_ext * mod$f_x))
    mod$norm$m_y + mod$norm$s_y * y_norm
}
