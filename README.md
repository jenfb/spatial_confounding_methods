---
title: 'Spatial confounding for individual-level exposures: sample code to fit the E-PS approach and other spatial models'
author: "Jennifer F. Bobb"
date: "2021-06-10"
output:
  html_document:
    keep_md: true
    theme: cerulean
    code_folding: hide
---



In this file, we provide sample code to fit the different spatial models considered in the paper 'Accounting for spatial confounding in epidemiological studies with individual-level exposures: an exposure penalized spline approach'. We run the models on a single generated dataset for one of the simulation scenarios considered in the paper.


```r
library(dplyr)
library(tidyr)
library(magrittr)

library(spBayes)
library(spNNGP)
library(mgcv)

source("functions/generate_spatial_data.R")
source("functions/extract_ests_gam.R")
```

## Models to compare

We consider the same set of models as from the simulation study, with the following exceptions (to make this sample code faster to run):

- we exclude the spline-based methods with a basis dimension of 1000
- for NNGP we use 5 rather than 10 nearest neighbors and run on 4 cores


```r
n_mcmc <- 10000
model_info <- bind_rows(
    tibble(method = "NS", inference = "freq"),
    tibble(method = "fix_df", inference = "freq", n_dim = c(5, 10, 25, 50, 100, 250, 500)),
    tibble(method = "PS", inference = "freq", n_dim = 500),
    tibble(method = "EPS", inference = "freq", n_dim = 500),
    tibble(method = "NNGP", inference = "bayes", n_mcmc = n_mcmc, n_neighbors = 5)
)
```

## Generate dataset from a particular simulation scenario


```r
## Parameters for all scenarios in simulation
#n_obs <- 2500
n_obs <- 1000
mat_corfun <- mat_corfun_matern
phi_seq <- c(0, 0.04, 0.15, 0.6)
params_all <- expand_grid(
    n_obs = n_obs,
    phi_c = phi_seq, phi_u = phi_seq,
    sigsq_x_true = c(0, round(1/18, 3), 1/2), 
    sigsq_y_true = c(3)^2,
    intercept = 0,
    beta_true = 3, 
    gamma_true = 1
)
params_all %<>% filter(!(sigsq_x_true == 0 & phi_u == 0))
## since if phi_u = 0 and sigsq_x = 0 then there is no variabilty in x
## that is not confounded (i.e., x is completely colinear with the confounder)
```


```r
## First, generate components used in data generating models across all scenarios
seed_dataset <- 10000
set.seed(seed_dataset)
df_comps <- gen_data_comps(
    n_obs = n_obs, phi_seq = phi_seq[phi_seq > 0],
    sigsq_x_true_seq = unique(params_all$sigsq_x_true),
    sigsq_y_true_seq = unique(params_all$sigsq_y_true),
    mat_corfun = mat_corfun,
    scale_proc = TRUE
)
#3.9 mins for n=2500 observations
#12 sec for n=1000 observations
##saveRDS(df_comps, file = "output/simulated_data_components.rds")
```


For illustration, here we generate a single dataset from one set of parameter values. To save computation time for this illustration we generate a dataset of size 1000 (in the paper we used 2500). The parameter values we use are defined here:

```r
params <- params_all %>%
    filter(phi_c == 0.04, phi_u == 0.60, sigsq_x_true == 0.5)
    
df <- params %$% gen_data_from_comps(
    df_comps, intercept = intercept, beta_true = beta_true, 
    sigsq_y_true = sigsq_y_true, sigsq_x_true = sigsq_x_true, 
    gamma_true = gamma_true, phi_c = phi_c, phi_u = phi_u
)
```



## Fit models to the simulated dataset

The following code loops through all of the methods to apply to the simulated dataset.


```r
res <- NULL
full_time0 <- Sys.time()
for(i in 1:nrow(model_info)) {
    message(i)
    res0 <- model_info[i, ]
    seed_iter <- seed_dataset
    
    if (res0$method == "NS") {
        time0 <- Sys.time()
        mod0 <- lm(y ~ x, data = df)
        s0 <- summary(mod0)
        cf <- "x"
        summ0 <- s0$coef[cf, ] %>%
            set_names(c("est", "se", "t", "p")) %>% as.list() %>%
            as_tibble() %>%
            mutate(param = "beta",
                   bic = BIC(mod0))
        time1 <- Sys.time()
        res_mod <- tibble(software = "lm", 
                          time0 = time0, time1 = time1) %>%
            crossing(summ0)
        
    } else if (res0$method == "EPS") {
        k <- res0$n_dim
        time0 <- Sys.time()
        mod0 <- gam(x ~ s(p0, p1, k = k), data = df)
        mod <- gam(y ~ x + s(p0, p1, k = k), data = df, sp = mod0$sp)
        s0 <- summary(mod)
        time1 <- Sys.time()
        res_mod <- extract_ests_gam(
            mod, s = s0, coef = "x", param = "beta") %>%
            mutate(software = "gam", time0 = time0, time1 = time1)
        
    } else if (res0$method == "PS") {
        k <- res0$n_dim
        time0 <- Sys.time()
        pspl <- gam(y ~ x + s(p0, p1, k = k), data = df)
        time1 <- Sys.time()
        res_mod <- extract_ests_gam(
            pspl, coef = "x", param = "beta") %>%
            mutate(software = "gam", time0 = time0, time1 = time1)
        
    } else if (res0$method == "fix_df") {
        DF <- res0$n_dim ##Note K = DF+1
        time0 <- Sys.time()
        mod <- gam(y ~ x + s(p0, p1, k = DF+1, fx = TRUE), data = df)
        s0 <- summary(mod)
        time1 <- Sys.time()
        res_mod <- extract_ests_gam(
            mod, s = s0, coef = "x", param = "beta") %>%
            mutate(software = "gam", time0 = time0, time1 = time1)
        
    } else if (res0$method == "NNGP") {
        ## set prior and tuning parmaters
        starting <- list("phi"=0.5, "sigma.sq"=1, "tau.sq"=1)
        tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
        priors <- list("beta.Flat",
                       "phi.Unif"=c(0.1, 30), "sigma.sq.IG"=c(0.01, 0.01),
                       "tau.sq.IG"=c(0.01, 0.01))
        cov.model <- "exponential"
        n.samples <- res0$n_mcmc
        #n.report <- 500
        burn.in <- 0.5*n.samples
        #burn.in <- 1
        
        ## fit model
        set.seed(seed_iter)
        time0 <- Sys.time()
        m.s <- df %$% spNNGP(y ~ x, 
                             coords = cbind(p0, p1), 
                             starting = starting, 
                             #method = "sequential", ## used in older version
                             method = "latent",
                             n.neighbors = res0$n_neighbors, tuning = tuning, 
                             priors = priors, cov.model = cov.model,
                             n.samples = n.samples, verbose=interactive(), 
                             #return.neighbors = FALSE, 
                             #n.omp.threads = 1
                             n.omp.threads = 4
        )
        time1 <- Sys.time()
        coef_info <- summary(m.s$p.beta.samples)$statistics["x", ]
        res_mod <- tibble(
            software = "spNNGP",
            est = coef_info["Mean"],
            se = coef_info["SD"],
            fixef_info = summary(m.s$p.beta.samples) %>% list(),
            vcomp_info = summary(m.s$p.theta.samples) %>% list(),
            time0 = time0, time1 = time1
        )
        
    } else {
        res_mod <- tibble(message = "method not available")
        
    }
    
    res_iter <- res0 %>%
        mutate(seed_iter = seed_iter) %>%
        bind_cols(res_mod)
     
    res %<>% bind_rows(res_iter)

}
full_time1 <- Sys.time()
#difftime(full_time1, full_time0) ## 1 min to run all the methods
res_all <- res
```

## Parameter estimates and computation time across models


```r
res %>% 
    mutate(seconds = as.numeric(difftime(time1, time0), units = "secs")) %>%
    transmute(method, `DF/dim.` = n_dim, beta_est = round(est, 2),
              `SE (or posterior SD)` = round(se, 3), EDF = round(`df_smooth`, 1),
              `runtime (sec)` = round(seconds, 2)
    )
```

```
## # A tibble: 11 x 6
##    method `DF/dim.` beta_est `SE (or posterior SD)`   EDF `runtime (sec)`
##    <chr>      <dbl>    <dbl>                  <dbl> <dbl>           <dbl>
##  1 NS            NA     3.37                  0.104    NA            0.01
##  2 fix_df         5     3.59                  0.113     5            0.08
##  3 fix_df        10     3.53                  0.118    10            0.07
##  4 fix_df        25     3.46                  0.12     25            0.1 
##  5 fix_df        50     3.37                  0.124    50            0.15
##  6 fix_df       100     3.09                  0.13    100            0.32
##  7 fix_df       250     2.87                  0.149   250            1.32
##  8 fix_df       500     2.69                  0.187   500            4.47
##  9 PS           500     3.51                  0.116    16           19.2 
## 10 EPS          500     3.16                  0.129   146           19.7 
## 11 NNGP          NA     3.4                   0.113    NA           17.4
```


