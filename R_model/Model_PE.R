#### Model position error to define best method ####
####################################################-

rm(list = ls())

## read in support functions (adjusted from Santon et al. 2023, 
# -> https://github.com/MSSanton/glmms_workflow)
source("./R_model/Linear modelling workflow_support functions.R") 

## read in helper functions (e.g. for rolling mean)
source("./help_functions.R") 

# increase memory size per worker
options(future.globals.maxSize = 1000 * 1024^2)

## load additional packages
library(dplyr)
library(sf)
library(ggeffects) # to predict gams (in fast way)
library(parallel)
#library(purrr) # kfold
#library(caret) # kfold
library(zoo) #rollM
library(geosphere) # distm()
library(ggspatial)
library(lubridate)

## global options
theme_set(theme_light()) # ggplot theme 
nsim <- 500 # number of simulations for distribution predictions
SAMP <- FALSE # use sample (nsamp, TRUE) or full data (FALSE)
nsamp <- 20000 # number of samples to use (subsample of whole data for faster modelling)

crsLL <- 4326 # coordinates in lon lat

## load Lookup table for models
# df.mod <- read.csv("./data/Lookup_model.csv")

## load data
dsn <- paste0("./data/cali/savedFiles/Data_cali_density.gpkg")
dat <- st_read(layer = 'density_est', dsn = dsn)
dat <-  dat %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

## load raster (to test animal data)
df.r <- st_read(dsn = "./data/cali/savedFiles/Data_cali_raster.gpkg", 
                layer = "density_raster")
df.r <- st_transform(df.r, crs = crsLL)

## check data structure, column names and vector classes. Change them as needed.
str(dat)

dat$r <- gsub('m', '', dat$r)
dat$r <- as.numeric(dat$r)

colnames(dat)[colnames(dat) == "Signal.max"] <- "maxSig"
colnames(dat)[colnames(dat) == "dens"] <- "cover"
colnames(dat)[colnames(dat) == "Station.Count"] <- "Sc"
colnames(dat)[colnames(dat) == "Antenna.Count"] <- "Ac"

## add tag height
dat$tagID[dat$Individual %in% c("TT090C", "TT241D")] <- 1.5
dat$tagID[dat$Individual %in% c("TT163C", "TT298D", "TT164D", "TT014D")] <- 1.0
dat$tagID[dat$Individual %in% c("TT240C", "TT090D")] <- 0.5
dat$tagID[dat$Individual %in% c("TT298C")] <- 2.0

## empty dfs for storing all results
df.coef <- NULL
df.pred <- NULL
df.sim <- NULL

#### 1) m1 to get r for meth = direct.ab ------------------------------------####
#------------------------------------------------------------------------------#
m <- "direct.ab"

for(s in c("maisC", "maisD")) {

  ## subset data
  df <- subset(dat, site == s & meth == m)
  
  ## convert columns to factors or ordered factors
  df <- df %>% mutate(across(c(tagID, meth), as.factor))
  df$r <- factor(df$r, ordered = T)
  
  ## remove NAs
  df <- df %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(r) & !is.na(cover))

  # ## use sample only (remove for final version)
  if(SAMP & nsamp <= nrow(df)) df <- df[sample(c(rep(TRUE, nsamp), rep(FALSE, nrow(df)-nsamp)), nrow(df) ,replace = F),]
  
  pdf(paste0("./output_model/output_m1_", s, "_", m, ".pdf"), height = 10, width = 10)
  
  print(ggplot(df) + geom_histogram(aes(x = PE)) + facet_wrap(~r) + ggtitle(paste(s, " - ", m)))
  
  ## model
  mod <- glmmTMB(PE ~ r + (1|tagID),
                 dispformula = ~ r,
                 family = lognormal(link = "log"),
                 # family = Gamma(link = "log"),
                 data = df)
  saveRDS(mod, paste0("./output_model/model_m1_", s, "_", m, ".RDS"))
  
  ## check residual plots (functions from Santon et al. 2023)
  # residual_plots(data = df, modelTMB = mod, response = "PE")
  # residual_plots_predictors(data = df, modelTMB = mod, predictors = "r")
  
  ## check dispersion parameters using raw (observed) and simulated values
  # disp.sim(data = df, modelTMB = mod, response = "PE", predictor = "r", n.sim = nsim)
  sd.sim(data = df, modelTMB = mod, response = "PE", predictor = "r", n.sim = nsim)
  mean.sim(data = df, modelTMB = mod, response = "PE", predictor = "r", n.sim = nsim)
  median.sim(data = df, modelTMB = mod, response = "PE", predictor = "r", n.sim = nsim)
  ppcheck_fun(data = df, modelTMB = mod, response = "PE", predictor = "r", n.sim = nsim)
  
  ## model coefficients and R2
  tmp.coef <- comp_int(modelTMB = mod,
           ci_range = 0.95,
           effects = "all",  
           component = "all")
  tmp.coef$model <- "m1"
  tmp.coef$site <- s
  tmp.coef$meth <- m
  
  df.coef <- bind_rows(df.coef, tmp.coef)
  
  # pseudo_r_squared
  1 - sum((mod$frame$PE - predict(mod, type = "response"))^2)/
    sum((mod$frame$PE - mean(mod$frame$PE))^2)

  newdat <- expand.grid("r" = factor(unique(df$r), ordered = T),
                        "tagID" = NA)
  
  # get model prediction (summarized)
  tmp.pred <- post_predictN(data = df, modelTMB = mod,
                           sims = TRUE, nsim = 4000,# save model simulations in mod.sims
                           newdat = newdat, component = "all",
                           DISP = T)
  tmp.pred$model <- "m1"
  tmp.pred$site <- s
  tmp.pred$meth <- m
  
  ## simulate raw data an get quantile (here: median)
  q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], 
                      probs = 0.5, nsim = 1000, SPEED = T)
  tmp.pred$q50 <- rowMeans(q50)
  tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
  tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
  tmp.pred$lwr25q50 <- apply(q50, 1, quantile, probs = 0.25)
  tmp.pred$upr75q50 <- apply(q50, 1, quantile, probs = 0.75)
  
  df.pred <- bind_rows(df.pred, tmp.pred)
  
  tmp.sim <- data.frame()
  for(i in 1:nrow(mod.sims[[1]])) {
    tmp.sim <- rbind(tmp.sim, data.frame(
      sim.m = mod.sims[[1]][i,],
      sim.sd = mod.sims[[2]][i,],
      sim.q50 = q50[i,],
      r = newdat$r[i]))
  }
  tmp.sim$model <- "m1"
  tmp.sim$site <- s
  tmp.sim$meth <- m
  
  df.sim <- bind_rows(df.sim, tmp.sim)
  
  ## final plot
  print(ggplot(tmp.sim) + 
    stat_halfeye(aes(x = r, y = sim.m,
                     linewidth = after_stat(.width)), # needed for linewidth
                 .width = c(0.5, 0.95),
                 color = "black",
                 fatten_point = 5) +
    scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
    # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
    ggtitle(paste0(s, " - ", m), 
            subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
    ylab("est. mean position error") +
    theme_light(base_size = 14) +
    theme(legend.position = "none"))
  
  print(ggplot(tmp.sim) + 
          stat_halfeye(aes(x = r, y = sim.q50,
                           linewidth = after_stat(.width)), # needed for linewidth
                       .width = c(0.5, 0.95),
                       color = "black",
                       fatten_point = 5) +
          scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
          # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
          ggtitle(paste(s, " - ", m), 
                  subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
          ylab("est. median position error") +
          theme_light(base_size = 14) +
          theme(legend.position = "none"))
  
  
  ## pairwise comparison
  pairwise_comparisons(data = df, modelTMB = mod, predictors = "r",
                       component = "cond", dispScale = "response", contrasts = "all")

  df %>% group_by(r) %>% 
    summarise(mean = mean(PE), median = median(PE), SD = sd(PE), n = n())
  
  dev.off()
  
  
}
write.csv(df.pred, "./output_model/model-predictions_m1.csv", row.names = F)
write.csv(df.sim, "./output_model/model-simulations_m1.csv", row.names = F)
write.csv(df.coef, "./output_model/model-coefficients_m1.csv", row.names = F)


## define selected values for r (based on models m1)
dat <- dat[!(dat$site == "maisC" & dat$meth == "direct.ab" & dat$r != 800),]
dat <- dat[!(dat$site == "maisD" & dat$meth == "direct.ab" & dat$r != 900),]

#### 2) m2 to get best meth ------------------------------------------------####
#------------------------------------------------------------------------------#

for(s in c("maisC", "maisD")) {
 
  df <- subset(dat, site == s)

  ## convert columns to factors or ordered factors
  df <- df %>% mutate(across(c(tagID, meth), as.factor))

  ## remove NAs
  df <- df %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(meth))
  
  write.csv(df, paste0("./output_model/data_m2_", s, ".csv"), row.names = F)
  
  ## only sample timestamps present in all methods
  df <- df %>% group_by(X_time, tagID) %>% # and same tagIDs???
    mutate(Nmeth = length(unique(meth)),
           nP = n()) %>%
    ungroup()

  df <- df[df$Nmeth == nlevels(df$meth),]

  # ## use sample only (remove for final version)
  if(SAMP & nsamp <= nrow(df)) df <- df[sample(c(rep(TRUE, nsamp), rep(FALSE, nrow(df)-nsamp)), nrow(df) ,replace = F),]

  pdf(paste0("./output_model/output_m2_", s, "_all.pdf"), height = 10, width = 10)

  print(ggplot(df) + geom_histogram(aes(x = PE)) + facet_wrap(~meth) + ggtitle(paste(s, " - all")))

  ## model
  mod <- glmmTMB(PE ~ meth + (1|tagID),
                 family = lognormal(link = "log"),
                 dispformula = ~ meth,
                 data = df)


  saveRDS(mod, paste0("./output_model/model_m2_", s, "_all.RDS"))

  ## check residual plots (functions from Santon et al. 2023)
  # residual_plots(data = df, modelTMB = mod, response = "PE")
  # residual_plots_predictors(data = df, modelTMB = mod, predictors = "Sc")

  ## check dispersion parameters using raw (observed) and simulated values
  # disp.sim(data = df, modelTMB = mod, response = "PE", predictor = "meth", n.sim = nsim)
  sd.sim(data = df, modelTMB = mod, response = "PE", predictor = "meth", n.sim = nsim)
  mean.sim(data = df, modelTMB = mod, response = "PE", predictor = "meth", n.sim = nsim)
  median.sim(data = df, modelTMB = mod, response = "PE", predictor = "meth", n.sim = nsim)
  ppcheck_fun(data = df, modelTMB = mod, response = "PE", predictor = "meth", n.sim = nsim)

  ## model coefficients and R2
  tmp.coef <- comp_int(modelTMB = mod,
                       ci_range = 0.95,
                       effects = "all",
                       component = "all")
  tmp.coef$model <- "m2"
  tmp.coef$site <- s

  df.coef <- bind_rows(df.coef, tmp.coef)

  # pseudo_r_squared
  1 - sum((mod$frame$PE - predict(mod, type = "response"))^2)/
    sum((mod$frame$PE - mean(mod$frame$PE))^2)

  newdat <- expand.grid("meth" = levels(df$meth),
                        "tagID" = NA)

  # get model prediction (summarized)
  tmp.pred <- post_predictN(data = df, modelTMB = mod,
                            sims = TRUE, nsim = 4000, # does not save model simulations in mod.sims
                            newdat = newdat,
                            DISP = T, 
                            component = "all")

  tmp.pred$model <- "m2"
  tmp.pred$site <- s


  ## get size of CI for plotting
  tmp.pred$diffCI <- tmp.pred$upr - tmp.pred$lwr

  ## simulate raw data an get quantile (here: median)
  q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], 
                      probs = 0.5, nsim = 1000, SPEED = T)
  tmp.pred$q50 <- rowMeans(q50)
  tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
  tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
  tmp.pred$lwr25q50 <- apply(q50, 1, quantile, probs = 0.25)
  tmp.pred$upr75q50 <- apply(q50, 1, quantile, probs = 0.75)
  
  df.pred <- bind_rows(df.pred, tmp.pred)
  
  ## mod.sim as long formate
  tmp.sim <- data.frame()
  for(i in 1:nrow(mod.sims[[1]])) {
    tmp.sim <- rbind(tmp.sim, data.frame(sim.m = mod.sims[[1]][i,],
                                         sim.sd = mod.sims[[2]][i,],
                                         sim.q50 = q50[i,],
                                         meth = newdat$meth[i]))
  }
  tmp.sim$model <- "m2"
  tmp.sim$site <- s

  df.sim <- bind_rows(df.sim, tmp.sim)

  ## final plot
  print(ggplot(tmp.sim) +
          stat_halfeye(aes(x = meth, y = sim.m,
                           linewidth = after_stat(.width)), # needed for linewidth
                       .width = c(0.5, 0.95),
                       color = "black",
                       fatten_point = 3) +
          scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
          # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
          ggtitle(paste(s, " - all"),
                  subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
          xlab("method") +
          ylab("est. mean position error") +
          theme_light(base_size = 14) +
          theme(legend.position = "none"))
  
  print(ggplot(tmp.sim) +
          stat_halfeye(aes(x = meth, y = sim.q50,
                           linewidth = after_stat(.width)), # needed for linewidth
                       .width = c(0.5, 0.95),
                       color = "black",
                       fatten_point = 3) +
          scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
          # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
          ggtitle(paste(s, " - all"),
                  subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
          xlab("method") +
          ylab("est. median position error") +
          theme_light(base_size = 14) +
          theme(legend.position = "none"))

  df %>% group_by(meth) %>%
    summarise(mean = mean(PE), median = median(PE), SD = sd(PE), n = n())
  dev.off()
  
}

write.csv(df.pred, "./output_model/model-predictions_m2.csv", row.names = F)
write.csv(df.sim, "./output_model/model-simulations_m2.csv", row.names = F)
write.csv(df.coef, "./output_model/model-coefficients_m2.csv", row.names = F)


#### 3) m3 to predict position error (pPE) ---------------------------------####
#------------------------------------------------------------------------------#
## exclude track for testing
df.test <- dat[dat$date %in% c("2023-10-13", "2023-09-13"),]
dat <- dat[!dat$date %in% c("2023-10-13", "2023-09-13"),]

df.coef <- NULL
df.pred <- NULL
df.pred2 <- NULL
df.sim <- NULL

for(s in c("maisC", "maisD")) {
  
  for(m in unique(dat$meth[dat$site == s])) {
    
    ## subset data
    df <- subset(dat, site == s & meth == m)
    
    ## test data
    df.t <- subset(df.test, site == s &  meth == m)
    
    ## convert columns to factors or ordered factors
    df <- df %>% mutate(across(c(tagID, meth), as.factor))
    
    ## get numeric predictors per method
    if(m %in% c("direct.ab", "omni.ab")) pred <- c("Sc", "cover", "Ac", "maxSig", "Weight")
    if(m %in% c("direct.an", "omni.ml")) pred <- c("Sc", "Ac", "cover", "maxSig")
    
    ## z-transform numeric variables
    df <- cbind(df, z_transform(data = df, predictors = pred))
    
    ## remove NAs
    df <- df %>% filter(if_all(all_of(c(pred, "PE", "tagID")), ~ !is.na(.)))
    df.t <- df.t %>% filter(if_all(all_of(c(pred, "PE", "tagID")), ~ !is.na(.)))
    
    pdf(paste0("./output_model/output_m5_", s, "_", m, ".pdf"), height = 10, width = 10)
    
    print(ggplot(df) + geom_histogram(aes(x = PE)) +
            geom_histogram(aes(x = PE), data = df.t, fill = "lightgrey") + 
            ggtitle(paste0(s, " - ", m), subtitle = "model data = dark, test data = light"))
    
    ## model formula per method
    
    if(m == "direct.ab") {
      mod.form <- formula(PE ~ Sc_z*Ac_z*cover_z + maxSig_z*Weight_z + (1|tagID))
      disp.form <- formula(~ Sc_z*Ac_z)
    }
    if(m == "omni.ab") {
      mod.form <- formula(PE ~ Sc_z*cover_z + maxSig_z*Weight_z + (1|tagID))
      disp.form <- formula(~ Sc_z)
    }
    if(m %in% c("direct.an", "omni.ml")) {
      mod.form <- formula(PE ~ Sc_z*cover_z + maxSig_z + (1|tagID))
      disp.form <- formula(~ Sc_z)
    }
    
    
    ## model
    mod <- glmmTMB(mod.form,
                   family = lognormal(link = "log"),
                   dispformula = disp.form,
                   data = df)
    
    saveRDS(mod, paste0("./output_model/model_m3_", s, "_", m, ".RDS"))
    
    ## simulate raw data
    xx <- simulate(mod, nsim = 1000)
    # write.csv(xx, paste0("./output_model/simulations_m3_", s, "_", m, ".csv"))
    df$meanPE <- apply(xx, 1, mean)
    df$q50PE <- apply(xx, 1, median)
    df$q65PE <- apply(xx, 1, quantile, probs = 0.65)
    df$q75PE <- apply(xx, 1, quantile, probs = 0.75)
    df$sdPE <- apply(xx, 1, sd)
    write.csv(df, paste0("./output_model/data_m3_", s, "_", m, ".csv"), row.names = F)
    
    ## check residual plots (functions from Santon et al. 2023)
    # residual_plots(data = df, modelTMB = mod, response = "PE")
    residual_plots_predictors(data = df, modelTMB = mod, predictors = "Sc")
    
    ## check dispersion parameters using raw (observed) and simulated values
    sd.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    mean.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    median.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    ppcheck_fun(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    
    ## model coefficients and R2
    tmp.coef <- comp_int(modelTMB = mod,
                         ci_range = 0.95,
                         effects = "all",
                         component = "all")
    tmp.coef$model <- "m3"
    tmp.coef$site <- s
    tmp.coef$meth <- m
    
    df.coef <- bind_rows(df.coef, tmp.coef)
    
    # pseudo_r_squared
    1 - sum((mod$frame$PE - predict(mod, type = "response"))^2)/
      sum((mod$frame$PE - mean(mod$frame$PE))^2)
    
    newdat <- expand.grid("Sc" = unique(df$Sc),
                          "Ac" = unique(df$Ac),
                          "tagID" = NA)
    
    newdat$Sc_z <- (newdat$Sc-mean(df$Sc))/sd(df$Sc)
    newdat$Ac_z <- (newdat$Ac-mean(df$Ac))/sd(df$Ac)
    
    ## get mean maxSig_z for each AcSc in raw data
    for(i in 1:nrow(newdat)) {
      Ac <- newdat$Ac[i]
      Sc <- newdat$Sc[i]
      newdat$maxSig_z[i] <- mean(df$maxSig_z[df$Ac == Ac & df$Sc == Sc])
      newdat$Weight_z[i] <- mean(df$Weight_z[df$Ac == Ac & df$Sc == Sc])
      newdat$cover_z[i] <- mean(df$cover_z[df$Ac == Ac & df$Sc == Sc])
    }
    
    newdat <- newdat[!is.nan(newdat$maxSig_z),] # remove NaNs due to non-present combinations
    
    ## predict for newdat
    tmp.pred <- post_predictN(data = df, modelTMB = mod,
                              sims = TRUE, nsim = 4000, # does save model simulations in mod.sims
                              newdat = dplyr::select(newdat, -c(Sc, Ac)),
                              DISP = T,
                              component = "all")
    ## predict for original data
    tmp.pred2 <- post_predictN(data = df, modelTMB = mod,
                               sims = F, nsim = 4000, # does not save model simulations in mod.sims
                               newdat = dplyr::select(df, -c(Sc, Ac, Weight, maxSig, cover)),
                               DISP = F,
                               component = "all")
    
    tmp.pred$model <- "m3"
    tmp.pred$site <- s
    tmp.pred$meth <- m
    
    tmp.pred2$model <- "m3"
    tmp.pred2$site <- s
    tmp.pred2$meth <- m
    
    ## simulate raw data an get quantile (here: median)
    # this is only needed if you want to predict something else than mean
    q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, 
                        nsim = 1000, SPEED = T)
    tmp.pred$q50 <- rowMeans(q50)
    tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
    tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
    tmp.pred$lwr25q50 <- apply(q50, 1, quantile, probs = 0.25)
    tmp.pred$upr75q50 <- apply(q50, 1, quantile, probs = 0.75)
    
    ## get size of CI for plotting
    tmp.pred$diffCI <- tmp.pred$upr - tmp.pred$lwr
    tmp.pred$diffCIq50 <- tmp.pred$uprq50 - tmp.pred$lwrq50
    
    df.pred <- bind_rows(df.pred, tmp.pred)
    df.pred2 <- bind_rows(df.pred2, tmp.pred2)
    
    ## mod.sim as long format
    tmp.sim <- data.frame()
    for(i in 1:nrow(mod.sims[[1]])) {
      tmp.sim <- rbind(tmp.sim, data.frame(sim.m = mod.sims[[1]][i,],
                                           sim.sd = mod.sims[[2]][i,],
                                           sim.q50 = q50[i,],
                                           Sc = newdat$Sc[i],
                                           Ac = newdat$Ac[i],
                                           maxSig_z = newdat$maxSig_z[i],
                                           Weight_z = newdat$Weight_z[i],
                                           cover_z = newdat$cover_z[i]))
    }
    tmp.sim$model <- "m3"
    tmp.sim$site <- s
    tmp.sim$meth <- m
    
    df.sim <- bind_rows(df.sim, tmp.sim)
    
    ## final plots
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = pred_mod, size = diffCI)) +
            scale_color_viridis_c("mean est. PE" ,option = "rocket", limits = c(0, 200), na.value = "#FAEBDDFF") +
            scale_size_continuous("range CI", range = c(2, 7)) +
            ggtitle(paste0(s, " - ", m)) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = q50, size = diffCIq50)) +
            scale_color_viridis_c("median est. PE" ,option = "rocket", limits = c(0, 200), na.value = "#FAEBDDFF") +
            scale_size_continuous("range CI", range = c(2, 7)) +
            ggtitle(paste0(s, " - ", m)) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = maxSig), size = 4) +
            scale_color_viridis_c(option = "rocket") +
            ggtitle(paste0(s, " - ", m)) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = Weight), size = 4) +
            scale_color_viridis_c(option = "rocket") +
            ggtitle(paste0(s, " - ", m)) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = Weight/Sc), size = 4) +
            scale_color_viridis_c(option = "rocket") +
            ggtitle(paste0(s, " - ", m)) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = cover), size = 4) +
            scale_color_viridis_c(option = "rocket") +
            ggtitle(paste0(s, " - ", m)) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.pred2) + 
            geom_point(aes(x = PE, y = pred_mod, color = Ac/Sc), pch = 1, size = 2, alpha = 0.4) + 
            #xlim(0, 500) + 
            geom_abline(slope = 1, lty = "dashed") + 
            scale_color_viridis_c("Ac/Sc") +
            facet_wrap(~Sc, scales = "free"))
    
    if(m %in% c("direct.ab", "omni.ab")) {
      print(ggplot(df) + geom_point(aes(x = Weight, y = maxSig, color = cover), alpha = 0.4) +
              scale_color_viridis_c(option = "rocket", direction = -1)) }
    
    if(m %in% c("direct.an", "omni.ml")) {
      print(ggplot(df) + geom_point(aes(x = cover, y = maxSig, color = Sc), alpha = 0.4) +
              scale_color_viridis_c(option = "rocket", direction = -1)) }
    
    df %>% group_by(Ac, Sc) %>%
      summarise(mean = mean(PE), median = median(PE), SD = sd(PE), n = n())
    
    ## predict for test testtrack
    ## scale numeric values
    newdat <- data.frame(tagID = rep(NA, nrow(df.t)))
    newdat$Sc_z <- (df.t$Sc-mean(df$Sc)) / sd(df$Sc)
    newdat$Ac_z <- (df.t$Ac-mean(df$Ac)) / sd(df$Ac)
    newdat$maxSig_z <- (df.t$maxSig-mean(df$maxSig)) / sd(df$maxSig)
    newdat$cover_z <- (df.t$cover-mean(df$cover)) / sd(df$cover)
    newdat$Weight_z <- (df.t$Weight-mean(df$Weight)) / sd(df$Weight)
    
    tmp.pred <- post_predictN(data = df, modelTMB = mod,
                              sims = TRUE, nsim = 1000,
                              newdat = newdat,
                              DISP = T,
                              component = "all")
    
    
    ## simulate raw data an get quantile (here: median)
    q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, 
                        nsim = 100, SPEED = T)
    tmp.pred$q50 <- rowMeans(q50)
    tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
    tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
    tmp.pred$lwr25q50 <- apply(q50, 1, quantile, probs = 0.25)
    tmp.pred$upr75q50 <- apply(q50, 1, quantile, probs = 0.75)
    
    ## get size of CI for plotting
    tmp.pred$diffCI <- tmp.pred$upr - tmp.pred$lwr
    tmp.pred$diffCIq50 <- tmp.pred$uprq50 - tmp.pred$lwrq50
    
    df.t <- cbind(df.t, tmp.pred[,colnames(tmp.pred)[!colnames(tmp.pred) %in% colnames(df.t)]]) ## check which rows to delete
    df.t$diff <- df.t$PE - df.t$pred_mod
    df.t$diffq50 <- df.t$PE - df.t$q50
    
    write.csv(df.t, paste0("./output_model/data_test_m3_", s, "_", m, ".csv"), row.names = F)
    
    print(ggplot(df.t) + 
            geom_abline(slope = 1, lty = "dashed", color = "blue") +
            geom_point(aes(x = PE, y = pred_mod), alpha = 0.2, pch = 1) +
            geom_linerange(aes(x = PE, ymin = lwr, ymax = upr), lwd = 0.5, alpha = 0.1))

    print(ggplot(df.t) + 
            geom_abline(slope = 1, lty = "dashed", color = "blue") +
            geom_point(aes(x = PE, y = q50), alpha = 0.2, pch = 1) +
            geom_linerange(aes(x = PE, ymin = lwrq50, ymax = uprq50), lwd = 0.5, alpha = 0.1))
    
    print(ggplot(df.t) + 
            stat_halfeye(aes(x = Individual, y = diff)))
    print(ggplot(df.t) + 
            stat_halfeye(aes(x = Individual, y = diffq50)))
    
    print(ggplot(df.t) + 
            stat_halfeye(aes(x = Sc, y = diff)))
    print(ggplot(df.t) + 
            stat_halfeye(aes(x = Sc, y = diffq50)))
    
    dev.off()
  
  }
}

write.csv(df.pred, "./output_model/model-predictions_m3.csv", row.names = F)
write.csv(df.pred2, "./output_model/model-predictions2_m3.csv", row.names = F)
write.csv(df.sim, "./output_model/model-simulations_m3.csv", row.names = F)
write.csv(df.coef, "./output_model/model-coefficients_m3.csv", row.names = F)

#### 4) test m3 with animal data -------------------------------------------####
#------------------------------------------------------------------------------#
## animal data (Great Tit, European Robin) from maisC tested with direct.ab
## get model and raw data
mod <- readRDS(paste0("./output_model/model_m3_maisC_direct.ab.RDS"))
df  <- read.csv("./output_model/data_m3_maisC_direct.ab.csv")
df.t <- st_read(dsn = "./data/animal/Data_animal_raw_maisC.gpkg")
shp.bird <- st_read(dsn = "./data/animal/Data_animal_handheld_maisC.gpkg")

shp.bird$Individual <- shp.bird$Testtag

# merge with df.r to get station cover (former 'dens')
df.t <- st_join(df.t, df.r[df.r$type == "direct",])
colnames(df.t)[colnames(df.t) == "dens"] <- "cover"

## scale numeric values
newdat <- data.frame(tagID = rep(NA, nrow(df.t)))
newdat$Sc_z <- (df.t$Sc-mean(df$Sc)) / sd(df$Sc)
newdat$Ac_z <- (df.t$Ac-mean(df$Ac)) / sd(df$Ac)
newdat$maxSig_z <- (df.t$maxSig-mean(df$maxSig)) / sd(df$maxSig)
newdat$cover_z <- (df.t$cover-mean(df$cover)) / sd(df$cover)
newdat$Weight_z <- (df.t$Weight-mean(df$Weight)) / sd(df$Weight)

## predict position errors
tmp.pred <- post_predictN(data = df, modelTMB = mod,
                          sims = TRUE, nsim = 4000,
                          newdat = newdat,
                          DISP = T,
                          component = "all")


# ## simulate raw data an get quantile (here: median and q65)
# q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, nsim = 1000)
# tmp.pred$q50 <- rowMeans(q50)
# tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
# tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
# 
# q65 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.65, nsim = 1000)
# tmp.pred$q65 <- rowMeans(q65)
# tmp.pred$lwrq65 <- apply(q65, 1, quantile, probs = 0.025)
# tmp.pred$uprq65 <- apply(q65, 1, quantile, probs = 0.975)

## get size of CI for plotting
tmp.pred$diffCI <- tmp.pred$upr - tmp.pred$lwr
# tmp.pred$diffCIq50 <- tmp.pred$uprq50 - tmp.pred$lwrq50
# tmp.pred$diffCIq65 <- tmp.pred$uprq65 - tmp.pred$lwrq65

df.t$pPE <- tmp.pred$pred_mod
# df.t$PA50 <- tmp.pred$q50
# df.t$PA65 <- tmp.pred$q65

## add date
df.t$date <- as.Date(df.t$X_time)

## save data
st_write(df.t, dsn = "./output_model/model_animal_m3_maisC_direct.ab.gpkg", 
         layer = "pPE", append = F)

# df.t <- st_read(dsn = "./output_model/model_animal_m3_maisC_direct.ab.gpkg",
#                 layer = "pPE")

## plot data
dim <- st_bbox(df.t)

print(ggplot(df.t[!is.na(df.t$lon.true),]) + 
  annotation_map_tile(type = "osm") +
  geom_sf(aes(color = pPE), alpha = 0.5, size = 3, pch = 1) +
  coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  scale_color_viridis_c(option = "inferno", limits = c(0, NA), end = 0.9, direction = -1) +
  facet_wrap(~Individual) +
  theme_void(base_size = 15))