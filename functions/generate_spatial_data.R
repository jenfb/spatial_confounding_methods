mat_corfun_matern <- function(phi, D) {
    Sigma <- (1+sqrt(3)*D/phi)*exp(-sqrt(3)*D/phi)
    Sigma
}

gen_data_comps <- function(n_obs = 100, phi_seq, sigsq_x_true_seq, sigsq_y_true_seq, rand_points = TRUE, mat_corfun = mat_corfun_exp, scale_proc = FALSE) {

    if (rand_points) {
        #n_obs=200
        locs <- matrix(runif(2*n_obs),ncol=2)
    } else {
        ## regular grid of points
        npt <- ceiling(sqrt(n_obs))
        sq <- (1:npt)/(npt + 1)
        locs <- crossing(p0 = sq, p1 = sq)
    }
    df <- locs %>%
        set_colnames(c("p0", "p1")) %>%
        as_tibble()

    D <- as.matrix(fields::rdist(as.matrix(locs)))
    proc_c <- proc_i <- matrix(NA, n_obs, length(phi_seq))
    for (j in seq_along(phi_seq)) {
        Sigma <- mat_corfun(phi = phi_seq[j], D = D)
        proc <- t(mvtnorm::rmvnorm(2,sigma=Sigma,method="svd"))
        if (scale_proc) proc %<>% scale()
        proc_c[, j] <- proc[, 1]
        proc_i[, j] <- proc[, 2]
    }
    colnames(proc_c) <- paste0("proc_c_", phi_seq)
    colnames(proc_i) <- paste0("proc_i_", phi_seq)
    proc_c %<>% as_tibble()
    proc_i %<>% as_tibble()

    sigsq_seq <- c(sigsq_x_true_seq, sigsq_y_true_seq)
    names(sigsq_seq) <- c(paste0("sigsq_x_", sigsq_x_true_seq),
                          paste0("sigsq_y_", sigsq_y_true_seq))
    mat <- matrix(NA, nrow = n_obs, ncol = length(sigsq_seq))
    colnames(mat) <- c(paste0("err_x_", sigsq_x_true_seq),
                       paste0("err_y_", sigsq_y_true_seq))
    for (j in seq_along(sigsq_seq)) {
        mat[, j] <- rnorm(n = n_obs, 0, sqrt(sigsq_seq[j]))
    }
    mat %<>% as_tibble()

    df %<>%
        bind_cols(proc_c) %>%
        bind_cols(proc_i) %>%
        bind_cols(mat)

    df_comps <- df
    df_comps
}
gen_data_from_comps <- function(df_comps, intercept, beta_true, phi_c, phi_u, gamma_true, sigsq_x_true, sigsq_y_true) {

    if (phi_c == 0) {
        z_c <- 0
    } else {
        z_c <- with(df_comps, get(paste0("proc_c_", phi_c)))
    }

    if (phi_u == 0) {
        z_i <- 0
    } else {
        z_i <- with(df_comps, get(paste0("proc_i_", phi_u)))
    }

    ## when there is either no independent or no confounded component of x, make it so that the total variability in x is the same as when there is both an independent and a confounded component
    coef_gen_x <- 0.5
    coef_gen_x <- ifelse(phi_c == 0 | phi_u == 0, sqrt(coef_gen_x), coef_gen_x)
    mu_x <- coef_gen_x*(z_c + z_i)

    eps_x <- with(df_comps, get(paste0("err_x_", sigsq_x_true)))
    eps_y <- with(df_comps, get(paste0("err_y_", sigsq_y_true)))

    df <- df_comps %>%
        mutate(z_c = z_c,
               z_i = z_i,
               mu_x = mu_x,
               #eps_x = rnorm(n = nrow(.), 0, sqrt(sigsq_x_true)),
               eps_x = eps_x,
               x = mu_x + eps_x,
               mu_y = intercept + beta_true*x + gamma_true*z_c,
               eps_y = eps_y,
               y = mu_y + eps_y
        ) %>%
        select(-contains("proc_"))

    df
}

