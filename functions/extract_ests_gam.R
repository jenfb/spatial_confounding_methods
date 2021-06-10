extract_ests_gam <- function(mod, modname = NULL, coef = "x", param = "beta", s = NULL) {
    if (is.null(s)) {
        s <- summary(mod)
    }

    #sp <- ifelse(is.null(mod$full.sp), ifelse(is.null(mod$sp), NA, mod$sp), NA)
    sp <- ifelse(is.null(mod$full.sp), ifelse(is.null(mod$sp), NA, mod$sp), mod$full.sp)
    summ0 <- s$p.table[coef, ] %>%
        set_names(c("est", "se", "t", "p")) %>%
        as.list() %>% as_tibble() %>%
        mutate(param = param,
               #df_smooth = s$s.table %>% as_tibble %$% edf,
               df_smooth = s$edf,
               #se_rob = sqrt(sandwich::sandwich(mod)[coef, coef]),
               bic = BIC(mod),
               gcv = mod$gcv.ubre,
               sp = sp)
    if (!is.null(modname)) {
        summ0 %<>% mutate(method = modname)
        summ0 %>%
            select(method, param, df_smooth, everything())
    } else {
        summ0 %>%
            select(param, df_smooth, everything())
    }

}
