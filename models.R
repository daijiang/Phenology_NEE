library(tidyverse)

# for all datasets, we will want to log10 transfer pop_km2
dat_agg = mutate(dat_agg, pop_km2_log10 = log10(pop_km2 + 1))

dat_agg %>% group_by(source, grp, sp) %>% tally() %>% arrange(n)

dat_agg %>% 
  group_by(source, grp) %>% 
  summarise(nspp = n_distinct(sp), 
            nsite = n_distinct(nth_cell),
            n = n())

mean(dat_agg$n_recd >= 5)

# clim_source = "worldclim"
# clim_source = "cru"
clim_source = "delaware"
# c("worldclim", "delaware", "cru")
temp_variable = paste("annu_temp_ave_11_5", clim_source, sep = "_")

dat_agg = dplyr::select(dat_agg, source, grp, nth_cell, sp, doy_median, n_recd, 
                        pop = pop_km2_log10, 
                        temp = temp_variable,
                        lsdt = LSDT_ave_2011_summer)
dat_agg = mutate(dat_agg, temp2 = temp^2, pop2 = pop^2)
dat_agg = na.omit(dat_agg)

envi = read_csv("data_output/local_popu_pep.csv") %>% 
  mutate(nth_cell = paste0("euro_", nth_cell)) %>% 
  bind_rows(read_csv("data_output/local_popu_npn_neon.csv") %>% 
              mutate(nth_cell = paste0("usa_", nth_cell))) %>% 
  dplyr::select(nth_cell, starts_with("LST")) %>% 
  filter(nth_cell %in% dat_agg$nth_cell) %>% 
  unique()

dat_agg = left_join(dat_agg, envi, by = "nth_cell")

summary(lm(lsdt ~ pop, data = dat_agg))

# which temp variable to use? change the formula below
fm_temp2_pop2_inter = doy_median ~ temp2 + pop2 + temp * pop + (1 | sp) + (1 | nth_cell) + (0 + temp | sp) + (0 + pop | sp)
fm_temp2_pop2 = doy_median ~ temp2 + pop2 +  + temp + pop + (1 | sp) + (1 | nth_cell) + (0 + temp | sp) + (0 + pop | sp)
fm_temp2 = doy_median ~ temp2 + temp + (1 | sp) + (1 | nth_cell) + (0 + temp | sp)
fm_pop2 = doy_median ~ pop2 + pop + (1 | sp) + (1 | nth_cell) + (0 + pop | sp)

# src = "all"
# gp = "flower"
# gp = "leaf"
lm_per_sp = msd = corr_results = mod_perse = mod_results = vector("list", length = 2)
i = 1
for(src in c("all")){ # c("PEP", "NPN_NEON", "all")
  for(gp in c("flower", "leaf")){
    nn = paste(gp, src, sep = "_")
    names(lm_per_sp)[i] = names(msd)[i] = names(mod_results)[i] = names(mod_perse)[i] = names(corr_results)[i] = nn
    cat(nn, "\n")
    
    if(src == "all"){
      dat = filter(dat_agg, grp == gp, temp < 16) # per reviewer's request
    } else {
      dat = filter(dat_agg, source == src, grp == gp)
    }
    
    num_cols = which(names(dat) == "n_recd") : ncol(dat)
    raw_msdat = bind_rows(
      dplyr::select(dat, num_cols) %>% summarise_all(mean),
      dplyr::select(dat, num_cols) %>% summarise_all(sd)
    ) # save mean and sd of raw data
    dat[, num_cols] = scale(dat[, num_cols])
    msd[[i]] = raw_msdat
    corr_results[[i]] = cor(dat[, num_cols])
    dat = mutate(dat, temppop = temp * pop)
    
    lm_per_sp[[i]] = group_by(dat, sp) %>% 
      do(broom::tidy(lm(doy_median ~ temp * pop, data = .))) %>% 
      ungroup()
# lme4::lmer
    if(gp == 'leaf' & src == 'all') cor_re = T else cor_re = F 
    # flower has problem with corr re and to speed up flower analysis
    # main results based on model weights do not change
    if(cor_re){
      fm_temp_pop_inter = doy_median ~ temp * pop + (1 + temp + pop + temp : pop | sp) + (1 | nth_cell)
      fm_temp_pop = doy_median ~ temp + pop + (1 + temp + pop|sp) + (1 | nth_cell)
      fm_temp = doy_median ~ temp + (1 + temp |sp) + (1 | nth_cell)
      fm_pop = doy_median ~ pop + (1 + pop |sp) + (1 | nth_cell)
    } else {
      fm_temp_pop_inter = doy_median ~ temp * pop + (1 | sp) + (0 + temp |sp) + (0 + pop|sp) + (0 + temp:pop | sp) + (1 | nth_cell)
      fm_temp_pop = doy_median ~ temp + pop + (1 | sp) + (0 + temp |sp) + (0 + pop|sp) + (1 | nth_cell)
      fm_temp = doy_median ~ temp + (1 | sp) +(0 + temp |sp) + (1 | nth_cell)
      fm_pop = doy_median ~ pop + (1 | sp) + (0 + pop|sp) + (1 | nth_cell)
    }
    
    dat = mutate(dat, continent = str_sub(nth_cell, 1, 4))
    
    mod_temp_pop_inter <- lmerTest::lmer(fm_temp_pop_inter, data = dat, REML = F)
    # mod_temp_pop_inter2 <- lmerTest::lmer(doy_median ~ temp * pop * continent + (1 | sp) + (0 + temp |sp) + (0 + pop|sp) + (0 + temp:pop | sp) + (1 | nth_cell), data = dat, REML = F)
    # summary(mod_temp_pop_inter)
    # anova(mod_temp_pop_inter2, mod_temp_pop_inter)
    
    mod_temp_pop <- lmerTest::lmer(fm_temp_pop, data = dat, REML = F)
    
    mod_pop <- lmerTest::lmer(fm_pop, data = dat, REML = F)
    
    mod_temp <- lmerTest::lmer(fm_temp, data = dat, REML = F)

    mod_perse[[i]] = list(mod_temp_pop_inter = mod_temp_pop_inter,
                          mod_temp_pop = mod_temp_pop,
                          mod_temp = mod_temp, mod_pop = mod_pop)
    
    sic_df = AIC(mod_temp_pop_inter, mod_temp_pop, mod_pop, mod_temp)
    mod_results[[i]] = purrr::map_df(c("mod_temp_pop_inter", "mod_temp_pop", "mod_pop", "mod_temp"), mod_summary) %>%
      spread(term, est) %>%
      rename(#temp_pop_inter = `temp:pop`,
             Intercept = `(Intercept)`) %>%
      dplyr::select(models, Intercept, temp, pop, `temp:pop`, everything()) %>%
      left_join(tibble(models = c("mod_temp_pop_inter", "mod_temp_pop", "mod_pop", "mod_temp"),
                           aic = sic_df$AIC,
                           mod_weights = as.vector(Weights(sic_df))),
                by = "models")
    i = i + 1
  }
}

rr2::R2(mod_perse$flower_all$mod_temp_pop_inter)
rr2::R2(mod_perse$flower_all$mod_temp_pop_inter, mod_perse$flower_all$mod_temp_pop) # contribution of pop:temp
rr2::R2(mod_perse$flower_all$mod_temp_pop_inter, mod_perse$flower_all$mod_temp) # contribution of pop
rr2::R2(mod_perse$flower_all$mod_temp_pop_inter, mod_perse$flower_all$mod_pop) # partial r2 of temp
rr2::R2(mod_perse$leaf_all$mod_temp_pop_inter)
rr2::R2(mod_perse$leaf_all$mod_temp_pop_inter, mod_perse$leaf_all$mod_temp_pop) # partial r2 of pop:temp
rr2::R2(mod_perse$leaf_all$mod_temp_pop_inter, mod_perse$leaf_all$mod_temp) # partial r2 of pop
rr2::R2(mod_perse$leaf_all$mod_temp_pop_inter, mod_perse$leaf_all$mod_pop) # partial r2 (contribution) of temp

mod_summary_list = function(lmm_list){
  plyr::ldply(lmm_list, function(x) {
    class(x) = 'merMod'
    broom::tidy(x)
  })
}
# summary of models
bind_rows(
  mod_summary_list(mod_perse$flower_all) %>% mutate(group = 'flower_all'),
  mod_summary_list(mod_perse$leaf_all) %>% mutate(group = 'leaf_all'))

# urban heat island ----
UHI_results = vector("list", 2)
names(UHI_results) = c("flower", "leaf")
for(gp in c('flower', 'leaf')){
  dat = filter(na.omit(dat_agg), grp == gp) %>% 
    mutate(uhi = UHI_11_5_delaware)
  num_cols = which(names(dat) == "n_recd") : ncol(dat)
  dat[, num_cols] = scale(dat[, num_cols])
  
  mod_pop_flower = lmerTest::lmer(doy_median ~ pop + (1 | sp) + (0 + pop | sp) + (1 | nth_cell), dat, REML = F)
  mod_temp_flower = lmerTest::lmer(doy_median ~ temp + (1 | sp) + (0 + temp | sp) + (1 | nth_cell), dat, REML = F)
  mod_temp_pop_flower = lmerTest::lmer(doy_median ~ temp + pop + (1 | sp) + (0 + temp | sp) + (0 + pop | sp) + (1 | nth_cell), dat, REML = F)
  MuMIn::r.squaredGLMM(mod_pop_flower)
  MuMIn::r.squaredGLMM(mod_temp_flower)
  MuMIn::r.squaredGLMM(mod_temp_pop_flower)
  
  mod_temp_lsdt_flower = lmerTest::lmer(doy_median ~ temp + uhi + (1 | sp) + (0 + temp |sp) + 
                                          (0 + uhi|sp) + (1 | nth_cell), data = dat, REML = F)
  anova(mod_temp_lsdt_flower, mod_temp_flower)
  mod_temp_lsdt_inter_flower = lmerTest::lmer(doy_median ~ temp * uhi + (1 | sp) + (0 + temp |sp) + 
                                                (0 + uhi|sp) + (0 + temp:uhi | sp) + 
                                                (1 | nth_cell), data = dat, REML = F)
  anova(mod_temp_lsdt_inter_flower, mod_temp_lsdt_flower)
  mod_temp_lsdt_pop_flower = lmerTest::lmer(doy_median ~ temp + uhi + pop + (1 | sp) + (0 + temp |sp) + 
                                              (0 + uhi|sp) + (0 + pop|sp) + (1 | nth_cell), data = dat, REML = F)
  mod_temp_inter_lsdt_pop_flower = lmerTest::lmer(doy_median ~ temp * uhi + pop + (1 | sp) + (0 + temp |sp) + 
                                                    (0 + uhi|sp)  + (0 + temp:uhi | sp) + 
                                                    (0 + pop | sp) + (1 | nth_cell), data = dat, REML = F)
  anova(mod_temp_inter_lsdt_pop_flower, mod_temp_lsdt_inter_flower)
  mod_temp_lsdt_pop_inter_flower = lmerTest::lmer(doy_median ~ temp + pop + uhi + 
                                                    temp:uhi + temp : pop + (1 | sp) + 
                                                    (0 + temp |sp) + (0 + uhi | sp) + 
                                                    (0 + temp:uhi|sp) + (0 + pop|sp) +
                                                    (0 + temp:pop | sp) + (1 | nth_cell), data = dat, REML = F)
  # mod_temp_lsdt_pop_inter_flower0 = lmerTest::lmer(doy_median ~ temp + pop + uhi + 
  #                                                    temp : pop + (1 | sp) + 
  #                                                    (0 + temp |sp) + (0 + uhi | sp) + 
  #                                                    (0 + pop|sp) +
  #                                                    (0 + temp:pop | sp) + (1 | nth_cell), data = dat, REML = F)
  aic_flower = AIC(mod_temp_flower, mod_temp_lsdt_flower, mod_temp_lsdt_inter_flower, mod_temp_inter_lsdt_pop_flower,
                   # mod_temp_lsdt_pop_inter_flower0, 
                   mod_temp_lsdt_pop_inter_flower)
  # plyr::ldply(c(mod_temp_flower, mod_temp_lsdt_flower, mod_temp_lsdt_inter_flower, mod_temp_inter_lsdt_pop_flower, mod_temp_lsdt_pop_inter_flower), MuMIn::r.squaredGLMM) %>% round(4)
  xt = lapply(c(mod_temp_flower, mod_temp_lsdt_flower, mod_temp_lsdt_inter_flower, mod_temp_inter_lsdt_pop_flower,
                # mod_temp_lsdt_pop_inter_flower0,
                mod_temp_lsdt_pop_inter_flower), broom.mixed::tidy, effects = "fixed")
  
  UHI_results[[gp]] = tibble(Models = c('Temp.', 'Temp. + LSDT', 'Temp. × LSDT', 'Temp. × LSDT + Pop.', 
                                        # 'Temp. + LSDT + Temp. × Pop.', 
                                        'Temp. × LSDT + Temp. × Pop.')) %>% 
    mutate(model_weights = Weights(aic_flower)) %>% 
    bind_cols(plyr::ldply(xt, function(x){
      mutate(x, estimate = paste0(round(estimate, 3), 
                                  ifelse(p.value < 0.001, "***", 
                                         ifelse(p.value < 0.01, "**",
                                                ifelse(p.value < 0.05, "*", 
                                                       ifelse(p.value < 0.1, ".", "")))))) %>%
        dplyr::select(term, estimate) %>% 
        spread(term, estimate)
    })) %>% 
    rename(Intercept = `(Intercept)`,
           `Temp.` = temp,
           UHI = uhi,
           `Temp.:UHI` = `temp:uhi`,
           Pop. = pop,
           `Temp.:Pop.` = `temp:pop`) %>% 
    bind_cols(lapply(c(mod_temp_flower, mod_temp_lsdt_flower, mod_temp_lsdt_inter_flower, mod_temp_inter_lsdt_pop_flower,
                       # mod_temp_lsdt_pop_inter_flower0,
                       mod_temp_lsdt_pop_inter_flower), MuMIn::r.squaredGLMM) %>% 
                plyr::ldply())
}


# ecoregions ----
lm_per_sp = msd = corr_results = mod_perse = mod_results = vector("list", length = 2)
i = 1
for(src in c("all")){
  for(gp in c("flower", "leaf")){
    nn = paste(gp, src, sep = "_")
    names(lm_per_sp)[i] = names(msd)[i] = names(mod_results)[i] = 
      names(mod_perse)[i] = names(corr_results)[i] = nn
    cat(nn, "\n")
    
    if(src == "all"){
      dat = filter(dat_agg, grp == gp)
    } else {
      dat = filter(dat_agg, source == src, grp == gp)
    }
    
    num_cols = which(names(dat) == "n_recd") : (ncol(dat) - 1)
    raw_msdat = bind_rows(
      dplyr::select(dat, num_cols) %>% summarise_all(mean),
      dplyr::select(dat, num_cols) %>% summarise_all(sd)
    )
    dat[, num_cols] = scale(dat[, num_cols])
    msd[[i]] = raw_msdat
    corr_results[[i]] = cor(dat[, num_cols])
    dat = mutate(dat, temppop = temp * pop)
    
    lm_per_sp[[i]] = group_by(dat, sp) %>% 
      do(broom::tidy(lm(doy_median ~ temp * pop, data = .))) %>% 
      ungroup()
    # lme4::lmer
    if(gp == 'leaf' & src == 'all') cor_re = T else cor_re = F # flower has problem with corr re
    if(cor_re){
      fm_temp_pop_inter = doy_median ~ temp * pop + (1 + temp + pop + temp : pop | sp) + (1 | nth_cell) + (1|ecoregion)
      fm_temp_pop = doy_median ~ temp + pop + (1 + temp + pop|sp) + (1 | nth_cell) + (1|ecoregion)
      fm_temp = doy_median ~ temp + (1 + temp |sp) + (1 | nth_cell) + (1|ecoregion)
      fm_pop = doy_median ~ pop + (1 + pop |sp) + (1 | nth_cell) + (1|ecoregion)
    } else {
      fm_temp_pop_inter0 = doy_median ~ temp * pop + (1 | sp) + (0 + temp |sp) + (0 + pop|sp) + (0 + temp:pop | sp) + (1 | nth_cell) + (1|ecoregion)
      fm_temp_pop_inter = doy_median ~ temp * pop + (1 | sp) + (0 + temp |sp) + (0 + pop|sp) + (0 + temp:pop | sp) + (0 + temp:pop | sp:ecoregion) + (1 | nth_cell) + (1|ecoregion)
      fm_temp_pop = doy_median ~ temp + pop + (1 | sp) + (0 + temp |sp) + (0 + pop|sp) + (1 | nth_cell) + (1|ecoregion)
      fm_temp = doy_median ~ temp + (1 | sp) +(0 + temp |sp) + (1 | nth_cell) + (1|ecoregion)
      fm_pop = doy_median ~ pop + (1 | sp) + (0 + pop|sp) + (1 | nth_cell) + (1|ecoregion)
    }
    
    # mod_temp_pop_inter0 <- lmerTest::lmer(fm_temp_pop_inter0, data = dat, REML = F)
    # summary(mod_temp_pop_inter0)
    mod_temp_pop_inter <- lmerTest::lmer(fm_temp_pop_inter, data = dat, REML = F)
    summary(mod_temp_pop_inter)
    
    mod_temp_pop <- lmerTest::lmer(fm_temp_pop, data = dat, REML = F)
    
    mod_pop <- lmerTest::lmer(fm_pop, data = dat, REML = F)
    
    mod_temp <- lmerTest::lmer(fm_temp, data = dat, REML = F)
    
    mod_perse[[i]] = list(mod_temp_pop_inter = mod_temp_pop_inter,
                          mod_temp_pop = mod_temp_pop,
                          mod_temp = mod_temp, mod_pop = mod_pop)
    
    sic_df = AIC(mod_temp_pop_inter, mod_temp_pop, mod_pop, mod_temp)
    mod_results[[i]] = purrr::map_df(c("mod_temp_pop_inter", "mod_temp_pop", "mod_pop", "mod_temp"), mod_summary) %>%
      spread(term, est) %>%
      rename(#temp_pop_inter = `temp:pop`,
        Intercept = `(Intercept)`) %>%
      dplyr::select(models, Intercept, temp, pop, `temp:pop`, everything()) %>%
      left_join(data_frame(models = c("mod_temp_pop_inter", "mod_temp_pop", "mod_pop", "mod_temp"),
                           aic = sic_df$AIC,
                           mod_weights = as.vector(Weights(sic_df))),
                by = "models")
    i = i + 1
  }
}
