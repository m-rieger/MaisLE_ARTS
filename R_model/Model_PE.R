#### Model position error to define best method ####
####################################################-

## CHANGE direct.in to direct.an ##
###################################
###################################
###################################
###################################

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
df.r <- st_read(dsn = "./data/cali/savedFiles/Data_cali_raster.gpkg", layer = "raster_site")
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

  ## subset data based on values in tmp.mod
  df <- subset(dat, site == s & meth == m)
  
  ## convert columns to factors or ordered factors
  df <- df %>% mutate(across(c(tagID, meth), as.factor))
  df$r <- factor(df$r, ordered = T)
  
  ## remove NAs
  df <- df %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(r) & !is.na(cover))

  # ## use sample only (remove for final version)
  if(SAMP & nsamp <= nrow(df)) df <- df[sample(c(rep(TRUE, nsamp), rep(FALSE, nrow(df)-nsamp)), nrow(df) ,replace = F),]
  
  pdf(paste0("./output_model/output_m1_", tmp.mod$site[m], "_", tmp.mod$meth[m], ".pdf"), height = 10, width = 10)
  
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


## define selected values for r (based on models)
# r
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


#### 3) m3 to predict position error (pPE) maisC -----------------------------####
#------------------------------------------------------------------------------#
## exclude track for testing
df.test <- dat[dat$date %in% c("2023-10-13", "2023-09-13"),]
dat <- dat[!dat$date %in% c("2023-10-13", "2023-09-13"),]

df.coef <- NULL
df.pred <- NULL
df.pred2 <- NULL
df.sim <- NULL

for(s in c("maisC", "maisD")) {
  m <- "direct.ab"
  
  ## subset data based on values in tmp.mod
  df <- subset(dat, site == s & meth == m)
  
  ## test data
  df.t <- subset(df.test, site == s &  meth == m)
  
  ## convert columns to factors or ordered factors
  df <- df %>% mutate(across(c(tagID, meth), as.factor))
  
  ## z-transform numeric variables
  df <- cbind(df, z_transform(data = df, predictors = c("Sc", "cover", 
                                                        "Ac", "maxSig",
                                                        "Weight")))
  
  ## remove NAs
  df <- df %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(Weight) & !is.na(maxSig) & 
                        !is.na(Sc) & !is.na(Ac) & !is.na(cover))
  df.t <- df.t %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(Weight) & !is.na(maxSig) & 
                            !is.na(Sc) & !is.na(Ac) & !is.na(cover))
  
  pdf(paste0("./output_model/output_m5_", s, "_", m, ".pdf"), height = 10, width = 10)
  
  print(ggplot(df) + geom_histogram(aes(x = PE)) +
          geom_histogram(aes(x = PE), data = df.t, fill = "lightgrey") + ggtitle(paste0(s, " - ", m)))
  
  ## model
  mod <- glmmTMB(PE ~ Sc_z*Ac_z*cover_z + maxSig_z*Weight_z + (1|tagID),
                 family = lognormal(link = "log"),
                 dispformula = ~ Sc_z*Ac_z,
                 # family = Gamma(link = "log"),
                 data = df)
  
  saveRDS(mod, paste0("./output_model/model_m3_", s, "_", m, ".RDS"))
  
  ## simulate raw data
  xx <- simulate(mod, nsim = 1000)
  write.csv(xx, paste0("./output_model/simulations_m3_", s, "_", m, ".csv"))
  df$meanPE <- apply(xx, 1, mean)
  df$q50PE <- apply(xx, 1, median)
  df$q65PE <- apply(xx, 1, quantile, probs = 0.65)
  df$q75PE <- apply(xx, 1, quantile, probs = 0.75)
  df$sdPE <- apply(xx, 1, sd)
  write.csv(df, paste0("./output_model/data_m3_", s, "_", m, ".csv"), row.names = F)
  
  ## check residual plots (functions from Santon et al. 2023)
  # residual_plots(data = df, modelTMB = mod, response = "PE")
  # residual_plots_predictors(data = df, modelTMB = mod, predictors = "Sc")
  
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
  
  df.coef <- bind_rows(df.coef, tmp.coef)
  
  # pseudo_r_squared
  1 - sum((mod$frame$PE - predict(mod, type = "response"))^2)/
    sum((mod$frame$PE - mean(mod$frame$PE))^2)
  
  newdat <- expand.grid("Sc" = unique(df$Sc),
                        "Ac" = unique(df$Ac),
                        "tagID" = NA)
  
  newdat$Sc_z <- (newdat$Sc-mean(df$Sc))/sd(df$Sc)
  newdat$Ac_z <- (newdat$Ac-mean(df$Ac))/sd(df$Ac)
  
  newdat <- newdat[newdat$Sc <= newdat$Ac &
                     newdat$Sc >= newdat$Ac/4,]
  
  ## get mean maxSig_z for each AcSc in raw data
  for(i in 1:nrow(newdat)) {
    Ac <- newdat$Ac[i]
    Sc <- newdat$Sc[i]
    newdat$maxSig_z[i] <- mean(df$maxSig_z[df$Ac == Ac & df$Sc == Sc])
    newdat$Weight_z[i] <- mean(df$Weight_z[df$Ac == Ac & df$Sc == Sc])
    newdat$cover_z[i] <- mean(df$cover_z[df$Ac == Ac & df$Sc == Sc])
  }
  
  newdat <- newdat[!is.nan(newdat$maxSig_z),] # remove NaNs due to non-present combinations
  
  # get model prediction (summarized)
  tmp.pred <- post_predictN(data = df, modelTMB = mod,
                            sims = TRUE, nsim = 4000, # does not save model simulations in mod.sims
                            newdat = dplyr::select(newdat, -c(Sc, Ac)),
                            DISP = T,
                            component = "all")
  ####################
  tmp.pred2 <- post_predictN(data = df, modelTMB = mod,
                            sims = F, nsim = 4000, # does not save model simulations in mod.sims
                            newdat = dplyr::select(df, -c(Sc, Ac, Weight, maxSig, cover)),
                            DISP = F,
                            component = "all")
  ######################
  
  tmp.pred$model <- "m3"
  tmp.pred$site <- s
  
  tmp.pred2$model <- "m3"
  tmp.pred2$site <- s
  
  ## simulate raw data an get quantile (here: median)
  q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, nsim = 1000, SPEED = T)
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
  
  df.sim <- bind_rows(df.sim, tmp.sim)
  
  ## final plot
  print(ggplot(tmp.pred) +
          geom_point(aes(x = Sc, y = Ac, color = pred_mod, size = diffCI)) +
          scale_color_viridis_c("mean est. PE" ,option = "rocket", limits = c(0, 200), na.value = "#FAEBDDFF") +
          scale_size_continuous("range CI", range = c(2, 7)) +
          ggtitle(paste0(s, " - ", m)) +
          ylim(1, 30) +
          xlim(1, 10) +
          theme_light(base_size = 14))
  
  print(ggplot(tmp.pred) +
          geom_point(aes(x = Sc, y = Ac, color = q50, size = diffCIq50)) +
          scale_color_viridis_c("median est. PE" ,option = "rocket", limits = c(0, 200), na.value = "#FAEBDDFF") +
          scale_size_continuous("range CI", range = c(2, 7)) +
          ggtitle(paste0(s, " - ", m)) +
          ylim(1, 30) +
          xlim(1, 10) +
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
  
  # ggplot(df) + geom_point(aes(x = Weight, y = PE), alpha = 0.4)
  # ggplot(df) + geom_point(aes(x = maxSig, y = PE), alpha = 0.4)
  # ggplot(df) + geom_point(aes(x = cover, y = PE), alpha = 0.4)
  # ggplot(df) + geom_point(aes(x = Weight, y = maxSig, color = PE), alpha = 0.4) +
  #   scale_color_viridis_c(option = "rocket", direction = -1) # best display
  
  print(ggplot(tmp.pred2) + 
    geom_point(aes(x = PE, y = pred_mod, color = round(Ac/Sc, 0)), pch = 1, size = 2, alpha = 0.4) + 
    xlim(0, 500) + 
    geom_abline(slope = 1) + 
    scale_color_viridis_c("Ac/Sc") +
    facet_wrap(~Sc))
  
  print(ggplot(df) + geom_point(aes(x = Weight, y = maxSig, color = cover), alpha = 0.4) +
    scale_color_viridis_c(option = "rocket", direction = -1)) # best display
  # ggplot(df) + geom_point(aes(x = Weight, y = cover, color = maxSig), alpha = 0.4) +
  #   scale_color_viridis_c(option = "rocket")
  # ggplot(df) + geom_point(aes(x = cover, y = maxSig, color = Weight), alpha = 0.4) +
  #   scale_color_viridis_c(option = "rocket")
  
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
  q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, nsim = 100, SPEED = T)
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

write.csv(df.pred, paste0("./output_model/model-predictions_m3_", m, ".csv"), row.names = F)
write.csv(df.pred2, paste0("./output_model/model-predictions2_m3_", m, ".csv"), row.names = F)
write.csv(df.sim, paste0("./output_model/model-simulations_m3_", m, ".csv"), row.names = F)
write.csv(df.coef, paste0("./output_model/model-coefficients_m3_", m, ".csv"), row.names = F)

## test with animal data
#-------------------------------------------------------------------------------

## get model and raw data
mod <- readRDS(paste0("./output_model/model_m3_maisC_direct.ab.RDS"))
df <- read.csv("./output_model/data_m3_maisC_direct.ab.csv")

df.bird <- read.csv("../data/data_test/GTdata_full_bird.csv")
df.bird$date_start <- as.POSIXct(df.bird$date_start, tz = "UTC")
df.bird$date_start <- with_tz(df.bird$date_start, "CET")
df.bird$date_stop <- as.POSIXct(df.bird$date_stop, tz = "UTC")
df.bird$date_stop <- with_tz(df.bird$date_stop, "CET")
colnames(df.bird)[colnames(df.bird) %in% c("lon", "lat")] <- c("lon.true", "lat.true")

## summarize df.bird by ID and transform to sf
shp.bird <- df.bird %>% group_by(lon.true, lat.true, ID, Testtag) %>%
  summarize(date_start = mean(date_start), .groups = "drop")
shp.bird <- st_as_sf(shp.bird, coords = c("lat.true", "lon.true"), crs = crsLL)

## loop through individuals
fL <- list.files(path = "../data/data_test/", pattern = "antennabeams.gpkg")

#####################
#####################
## add raw data of animals to repository
#####################
#####################

for(f in fL) {
  dsn <- paste0("../data/data_test/", f)
  #dsn = "../data/data_test/GT13_R2_S1_Test.antennabeams.gpkg"
  df.t <- st_read(dsn = dsn)
  # df.t <- st_read(dsn = dsn, layer = "PA")
  
  df.t <- st_transform(df.t, crs = crsLL)
  
  # merge using sf objects
  df.t <- st_join(df.t, df.r)
  
  # get original positions from manual tracking
  df.t <- left_join(df.t, df.bird[, c("Testtag", "date_start", "lat.true", "lon.true", "ID")], 
                    by = c("Individual" = "Testtag", "X_time" = "date_start"))
  
  # using a raster
  # df.t$cover   <- extract(df.r2, df.t)
  
  ## change colnames
  colnames(df.t)[colnames(df.t) %in% c("Station.Count", "Antenna.Count", "Signal.max", "dens")] <- c("Sc", "Ac", "maxSig", "cover")
  
  ## scale numeric values
  newdat <- data.frame(tagID = rep(NA, nrow(df.t)))
  newdat$Sc_z <- (df.t$Sc-mean(df$Sc)) / sd(df$Sc)
  newdat$Ac_z <- (df.t$Ac-mean(df$Ac)) / sd(df$Ac)
  newdat$maxSig_z <- (df.t$maxSig-mean(df$maxSig)) / sd(df$maxSig)
  newdat$cover_z <- (df.t$cover-mean(df$cover)) / sd(df$cover)
  newdat$Weight_z <- (df.t$Weight-mean(df$Weight)) / sd(df$Weight)
  
  tmp.pred <- post_predictN(data = df, modelTMB = mod,
                            sims = TRUE, nsim = 4000,
                            newdat = newdat,
                            DISP = T,
                            component = "all")
  
  
  # ## simulate raw data an get quantile (here: median)
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
  
  df.t$PA <- tmp.pred$pred_mod
  # df.t$PA50 <- tmp.pred$q50
  # df.t$PA65 <- tmp.pred$q65
  
  ## save data
  st_write(df.t, dsn = dsn, layer = "PA", append = F)
  df.t <- st_read(dsn = dsn, layer = "PA")
  shp.tmp <- shp.bird[shp.bird$Testtag == unique(df.t$Individual) &
                        shp.bird$date_start >= min(df.t$X_time) &
                        shp.bird$date_start <= max(df.t$X_time),]
  
  dim <- st_bbox(df.t)
  
  print(ggplot(df.t[!is.na(df.t$lon.true),]) + 
    annotation_map_tile(type = "osm") +
    geom_sf(aes(color = PA), alpha = 0.5, size = 3, pch = 1) +
    coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
    annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
    scale_color_viridis_c(option = "rocket", limits = c(0, NA), na.value = "#FAEBDDFF") +
    ggtitle(f) +
    theme_void(base_size = 15))
  
  print(ggplot(df.t[!is.na(df.t$lon.true),]) + 
    annotation_map_tile(type = "osm") +
    geom_sf(aes(size = PA, color = X_time), alpha = 0.5, pch = 1) +
    geom_sf(aes(fill = date_start, color = date_start, shape = "23"), alpha = 0.5, size = 10, color = "black",
            data = shp.tmp) +
    scale_size_continuous("PA est. points",range = c(2, 10)) +
    coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
    annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
    scale_color_viridis_c("time", option = "inferno", end = 0.8, begin = 0.1,
                          limits = c(min(df.t$X_time, shp.tmp$date_start), max(df.t$X_time, shp.tmp$date_start))) +
    scale_fill_viridis_c("time", option = "inferno", end = 0.8, begin = 0.1,
                         limits = c(min(df.t$X_time, shp.tmp$date_start), max(df.t$X_time, shp.tmp$date_start))) +
    ggtitle(f) +
    theme_void(base_size = 15) +
      guides(shape = guide_legend(override.aes = list(color = "black"))) + 
      scale_shape_manual(name = "'ground truth' \npoints", values = c(23), labels = c(" ")))

  # ggplot(df.t[df.t$Ac > df.t$Sc,]) + 
  #   annotation_map_tile(type = "osm") +
  #   geom_sf(aes(color = PA), alpha = 0.5, size = 3, pch = 1) +
  #   scale_size_continuous(range = c(2, 10)) +
  #   coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
  #   annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  #   scale_color_viridis_c(option = "rocket", limits = c(0, NA), na.value = "#FAEBDDFF") +
  #   theme_void(base_size = 15)
  # 
  # ggplot(df.t[df.t$Ac > 1.5*df.t$Sc,]) + 
  #   annotation_map_tile(type = "osm") +
  #   geom_sf(aes(color = PA), alpha = 0.5, size = 3, pch = 1) +
  #   scale_size_continuous(range = c(2, 10)) +
  #   coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
  #   annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  #   scale_color_viridis_c(option = "rocket", limits = c(0, NA), na.value = "#FAEBDDFF") +
  #   theme_void(base_size = 15)
  # 
  # ggplot(df.t) + geom_point(aes(x = Ac, y = Sc), alpha = 0.01, size = 3)
  
  nrow(df.t)/nrow(df.t)
  nrow(df.t[df.t$Ac > df.t$Sc,])/nrow(df.t)
  nrow(df.t[df.t$Ac > 1.5*df.t$Sc,])/nrow(df.t)
  
}


#### 5) m5 to get position accuracy (PA) maisD -----------------------------####
#------------------------------------------------------------------------------#


#tmp.mod <- df.mod[df.mod$model == "m4" & df.mod$site == "maisD" & df.mod$meth == methD,]


for(m in c("direct.in", "omni.ml", "omni.ab")) {
  df.coef <- NULL
  df.pred <- NULL
  df.pred2 <- NULL
  df.sim <- NULL  
  
  for(s in c("maisC", "maisD")) {
  
  
    if(s == "maisD" & m %in% c("omni.ab", "omni.ml")) next
    ## subset data based on values in tmp.mod
    df <- subset(dat, 
                 site == s &
                   meth == m)
    
    ## test data
    df.t <- subset(df.test,
                   site == s &
                     meth == m)
    
    ## convert columns to factors or ordered factors
    df <- df %>% mutate(across(c(tagID, meth), as.factor))
    
    ## z-transform numeric variables
    df <- cbind(df, z_transform(data = df, predictors = c("Sc", "cover", 
                                                          "Ac", "maxSig",
                                                          "Weight")))
    
    ## remove NAs
    df <- df %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(maxSig) & 
                          !is.na(Sc) & !is.na(cover))
    df.t <- df.t %>% filter(!is.na(PE) & !is.na(tagID) & !is.na(maxSig) & 
                              !is.na(Sc) & !is.na(cover))
    
    # ## use sample only (remove for final version)
    if(SAMP & nsamp <= nrow(df)) df <- df[sample(c(rep(TRUE, nsamp), rep(FALSE, nrow(df)-nsamp)), nrow(df) ,replace = F),]

    pdf(paste0("./output_model/output_m5_", s, "_", m, ".pdf"), height = 10, width = 10)
    
    print(ggplot(df) + geom_histogram(aes(x = PE)) + ggtitle(paste0(s, " - ", m)))
    
    ## kfold cross validation
    
    ## mod formula
    mod.form <- list()
    #mod.form[["mfull"]]  <- formula(PE ~ Sc_z*Ac_z*cover_z*maxSig_z*Weight_z + (1|tagID)) # too much
    mod.form[["mfull1"]] <- formula(PE ~ Sc_z*cover_z*maxSig_z + (1|tagID))
    mod.form[["mfull2"]] <- formula(PE ~ Sc_z*cover_z + maxSig_z + (1|tagID))
    mod.form[["mfull3"]] <- formula(PE ~ Sc_z + cover_z*maxSig_z + (1|tagID))
    mod.form[["m11"]] <- formula(PE ~ Sc_z + cover_z + maxSig_z + (1|tagID))
    mod.form[["m21"]] <- formula(PE ~ Sc_z*cover_z + (1|tagID))
    mod.form[["m22"]] <- formula(PE ~ Sc_z*maxSig_z + (1|tagID))
    mod.form[["m23"]] <- formula(PE ~ cover_z*maxSig_z + (1|tagID))
    mod.form[["m6"]]  <- formula(PE ~ (1|tagID))
    
    if(m == "omni.ab") {
      mod.form <- list()
      mod.form[["mfull2"]] <- formula(PE ~ Sc_z*cover_z + maxSig_z*Weight_z + (1|tagID))
      mod.form[["mfull3"]] <- formula(PE ~ Sc_z*cover_z + maxSig_z + Weight_z + (1|tagID))
      mod.form[["m11"]] <- formula(PE ~ Sc_z + cover_z + maxSig_z + Weight_z + (1|tagID))
      mod.form[["m21"]] <- formula(PE ~ Sc_z*cover_z + maxSig_z + (1|tagID))
      mod.form[["m22"]] <- formula(PE ~ Sc_z*cover_z + Weight_z + (1|tagID))
      mod.form[["m23"]] <- formula(PE ~ Sc_z + maxSig_z + Weight_z + (1|tagID))
      mod.form[["m31"]] <- formula(PE ~ Sc_z*cover_z + (1|tagID))
      mod.form[["m32"]] <- formula(PE ~ Sc_z + Weight_z + (1|tagID))
      mod.form[["m33"]] <- formula(PE ~ Sc_z + maxSig_z + (1|tagID))

      mod.form[["m6"]]  <- formula(PE ~ (1|tagID))
      
    }
    # mod.form[["m31"]] <- formula(PE ~ Sc_z + cover_z + (1|tagID))
    # mod.form[["m32"]] <- formula(PE ~ Sc_z + maxSig_z + (1|tagID))
    # mod.form[["m33"]] <- formula(PE ~ cover_z + maxSig_z + (1|tagID))
    # mod.form[["m41"]]  <- formula(PE ~ Sc_z + (1|tagID))
    # mod.form[["m42"]]  <- formula(PE ~ cover_z + (1|tagID))
    # mod.form[["m43"]]  <- formula(PE ~ maxSig_z + (1|tagID))
    
    if(FOLD) {
      # Create k-folds
      set.seed(3171)  # For reproducibility
      folds <- createFolds(df$PE, k = k, list = TRUE)
      
      
      # Function to perform k-fold cross-validation
      k_fold_cv <- function(data, folds, mod.form, meth = "MAE") {
        results <- map_dbl(folds, function(test_indices) {
          # Split the data into training and testing sets
          train_data <- data[-test_indices, ]
          test_data <- data[test_indices, ]
          
          # Fit the model on the training data
          model <- glmmTMB(mod.form, data = train_data, 
                           dispformula = ~ Sc_z, family = lognormal(link = "log"))
          
          # Make predictions on the test data
          predictions <- predict(model, newdata = test_data, type = "response", re.form = NULL)
          
          # Calculate performance metric (e.g., RMSE)
          if(meth == "MAE") {
            mae <- mean(abs(test_data$PE - predictions))
            return(mae)          
          }
          
          if(meth == "RMSE") {
            rmse <- sqrt(mean((test_data$PE - predictions)^2))
            return(rmse)          
          }
          
        })
        
        return(results)
      }
      
      
      # Run k-fold cross-validation
      df.k <- data.frame(k = 1:k)
      
      for(i in names(mod.form)) df.k[, i] <- k_fold_cv(df, folds, mod.form[[i]])
      
      write.csv(df.k, paste0("./output_model/kfold_m5_", s, "_", m, ".csv"), row.names = F)
      colMeans(df.k)
      
      ## get differences and se diff
      k.sum <- data.frame(model = colnames(df.k)[colnames(df.k) != "k"],
                          RMSE = NA, RMSE_diff = NA, se_diff = NA, lwr = NA, upr = NA)
      
      for(i in k.sum$model) {
        k.sum$RMSE[k.sum$model == i] <- mean(df.k[, i])
        k.sum$RMSE_diff[k.sum$model == i] <- mean(df.k[, "mfull2"] - df.k[, i])
        k.sum$se_diff[k.sum$model == i] <- sd(df.k[, "mfull2"] - df.k[, i])/sqrt(nrow(df.k))
        k.sum$lwr[k.sum$model == i] <- quantile(df.k[, "mfull2"] - df.k[, i], probs = 0.025)
        k.sum$upr[k.sum$model == i] <- quantile(df.k[, "mfull2"] - df.k[, i], probs = 0.975)
        
      }
    }
    
    
    ## run best model
    ## model
    mod <- glmmTMB(mod.form[["mfull2"]],
                   family = lognormal(link = "log"),
                   dispformula = ~ Sc_z,
                   # family = Gamma(link = "log"),
                   data = df)
    
    
    saveRDS(mod, paste0("./output_model/model_m5_", s, "_", m, ".RDS"))
    
    ## simulate raw data
    xx <- simulate(mod, nsim = 1000)
    write.csv(xx, paste0("./output_model/simulations_m5_", s, "_", m, ".csv"))
    df$meanPE <- apply(xx, 1, mean)
    df$q50PE <- apply(xx, 1, median)
    df$q65PE <- apply(xx, 1, quantile, probs = 0.65)
    df$q75PE <- apply(xx, 1, quantile, probs = 0.75)
    df$sdPE <- apply(xx, 1, sd)
    write.csv(df, paste0("./output_model/data_m5_", s, "_", m, ".csv"), row.names = F)
    
    ## check dispersion parameters using raw (observed) and simulated values
    # disp.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    sd.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    mean.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    median.sim(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    ppcheck_fun(data = df, modelTMB = mod, response = "PE", n.sim = nsim)
    
    ## model coefficients and R2
    tmp.coef <- comp_int(modelTMB = mod,
                         ci_range = 0.95,
                         effects = "all",
                         component = "all")
    tmp.coef$model <- "m5"
    tmp.coef$site <- s
    
    df.coef <- bind_rows(df.coef, tmp.coef)
    
    # pseudo_r_squared
    1 - sum((mod$frame$PE - predict(mod, type = "response"))^2)/
      sum((mod$frame$PE - mean(mod$frame$PE))^2)
    
    r2(model = mod)
    
    newdat <- expand.grid("Sc" = unique(df$Sc),
                          "Ac" = unique(df$Ac),
                          "tagID" = NA,
                          "Weight_z" = NA)
    
    newdat$Sc_z <- (newdat$Sc-mean(df$Sc))/sd(df$Sc)
    newdat$Ac_z <- (newdat$Ac-mean(df$Ac))/sd(df$Ac)
    
    if(m == "direct.in") newdat <- newdat[2*newdat$Sc == newdat$Ac,]
    if(m != "direct.in") newdat <- newdat[newdat$Sc == newdat$Ac,]
    
    ## get mean maxSig_z for each AcSc in raw data
    for(i in 1:nrow(newdat)) {
      Ac <- newdat$Ac[i]
      Sc <- newdat$Sc[i]
      if(m == "omni.ab") newdat$Weight_z[i] <- mean(df$Weight_z[df$Ac == Ac & df$Sc == Sc])
      newdat$maxSig_z[i] <- mean(df$maxSig_z[df$Ac == Ac & df$Sc == Sc])
      newdat$cover_z[i] <- mean(df$cover_z[df$Ac == Ac & df$Sc == Sc])
    }
    
    newdat <- newdat[!is.nan(newdat$maxSig_z),] # remove NaNs due to non-present combinations
    
    
    tmp.pred <- post_predictN(data = df, modelTMB = mod,
                              sims = TRUE, nsim = 4000, # does not save model simulations in mod.sims
                              newdat = dplyr::select(newdat, -c(Sc, Ac)),
                              DISP = T,
                              component = "all")
    
    tmp.pred2 <- post_predictN(data = df, modelTMB = mod,
                              sims = F, nsim = 4000, # does not save model simulations in mod.sims
                              newdat = dplyr::select(df, -c(Sc, Ac, Weight, maxSig, cover)),
                              DISP = F,
                              component = "all")
    
    ## simulate raw data an get quantile (here: median)
    q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, nsim = 1000)
    tmp.pred$q50 <- rowMeans(q50)
    tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
    tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
    
    q65 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.65, nsim = 1000)
    tmp.pred$q65 <- rowMeans(q65)
    tmp.pred$lwrq65 <- apply(q65, 1, quantile, probs = 0.025)
    tmp.pred$uprq65 <- apply(q65, 1, quantile, probs = 0.975)
    
    ## get size of CI for plotting
    tmp.pred$diffCI <- tmp.pred$upr - tmp.pred$lwr
    tmp.pred$diffCIq50 <- tmp.pred$uprq50 - tmp.pred$lwrq50
    tmp.pred$diffCIq65 <- tmp.pred$uprq65 - tmp.pred$lwrq65
    
    tmp.pred$model <- "m5"
    tmp.pred$site <- s
    df.pred <- bind_rows(df.pred, tmp.pred)
    
    tmp.pred2$model <- "m5"
    tmp.pred2$site <- s
    df.pred2 <- bind_rows(df.pred2, tmp.pred2)
    
    ## mod.sim as long formate
    tmp.sim <- data.frame()
    for(i in 1:nrow(mod.sims[[1]])) {
      tmp.sim <- rbind(tmp.sim, data.frame(sim.m = mod.sims[[1]][i,],
                                           sim.sd = mod.sims[[2]][i,],
                                           sim.q50 = q50[i,],
                                           sim.q65 = q65[i,],
                                           Sc = newdat$Sc[i],
                                           Ac = newdat$Ac[i],
                                           maxSig_z = newdat$maxSig_z[i],
                                           #Weight_z = newdat$Weight_z[i],
                                           cover_z = newdat$cover_z[i]))
    }
    tmp.sim$model <- "m5"
    tmp.sim$site <- s
    
    df.sim <- bind_rows(df.sim, tmp.sim)
    
    ## final plot
    print(ggplot(tmp.pred) +
            geom_point(aes(x = Sc, y = Ac, color = median, size = diffCI)) +
            scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
            scale_size_continuous("range CI", range = c(2, 7)) +
            ggtitle(paste0(s, " - ", m)) +
            ylim(1, 20) +
            xlim(1, 10) +
            theme_light(base_size = 14))
    
    print(ggplot(tmp.sim) +
            stat_halfeye(aes(x = Sc, y = sim.m,
                             linewidth = after_stat(.width)), # needed for linewidth
                         .width = c(0.5, 0.95),
                         color = "black",
                         fatten_point = 3) +
            scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
            # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
            ggtitle(paste0(s, " - ", m),
                    subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
            xlab("method") +
            ylab("est. mean position error") +
            theme_light(base_size = 14) +
            theme(legend.position = "none"))
    
    print(ggplot(tmp.sim) +
            stat_halfeye(aes(x = Sc, y = sim.q50,
                             linewidth = after_stat(.width)), # needed for linewidth
                         .width = c(0.5, 0.95),
                         color = "black",
                         fatten_point = 3) +
            scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
            # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
            ggtitle(paste0(s, " - ", m),
                    subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
            xlab("method") +
            ylab("est. median position error") +
            theme_light(base_size = 14) +
            theme(legend.position = "none"))
    
    print(ggplot(tmp.sim) +
            stat_halfeye(aes(x = Sc, y = sim.q65,
                             linewidth = after_stat(.width)), # needed for linewidth
                         .width = c(0.5, 0.95),
                         color = "black",
                         fatten_point = 3) +
            scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
            # scale_linewidth_manual(values = c(`0.5` = 0.8, `0.95` = 1.5)) +
            ggtitle(paste0(s, " - ", m),
                    subtitle = "mean with 50% (thick) and 95% (thin) CI ") +
            xlab("method") +
            ylab("est. q65 position error") +
            theme_light(base_size = 14) +
            theme(legend.position = "none"))
    
    print(ggplot(tmp.pred2) + 
      geom_point(aes(x = PE, y = pred_mod, color = round(Ac/Sc, 0)), pch = 1, size = 2, alpha = 0.4) + 
      xlim(0, 500) + 
      geom_abline(slope = 1) + 
      scale_color_viridis_c("Ac/Sc") +
      facet_wrap(~as.factor(round(Sc, 0))))
    
    
    ## predict for test testtrack
    ## scale numeric values
    newdat <- data.frame(tagID = rep(NA, nrow(df.t)))
    newdat$Sc_z <- (df.t$Sc-mean(df$Sc)) / sd(df$Sc)
    newdat$maxSig_z <- (df.t$maxSig-mean(df$maxSig)) / sd(df$maxSig)
    newdat$cover_z <- (df.t$cover-mean(df$cover)) / sd(df$cover)
    if(m == "omni.ab") newdat$Weight_z <- (df.t$Weight-mean(df$Weight)) / sd(df$Weight)
    
    tmp.pred <- post_predictN(data = df, modelTMB = mod,
                              sims = TRUE, nsim = 1000,
                              newdat = newdat,
                              DISP = T,
                              component = "all")
    
    
    ## simulate raw data an get quantile (here: median)
    q50 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.5, nsim = 100, SPEED = T)
    tmp.pred$q50 <- rowMeans(q50)
    tmp.pred$lwrq50 <- apply(q50, 1, quantile, probs = 0.025)
    tmp.pred$uprq50 <- apply(q50, 1, quantile, probs = 0.975)
    tmp.pred$lwr25q50 <- apply(q50, 1, quantile, probs = 0.25)
    tmp.pred$upr75q50 <- apply(q50, 1, quantile, probs = 0.75)
    
    q65 <- quant.rlnorm(m = mod.sims[[1]], sd = mod.sims[[2]], probs = 0.65, nsim = 100, SPEED = T)
    tmp.pred$q65 <- rowMeans(q65)
    tmp.pred$lwrq65 <- apply(q65, 1, quantile, probs = 0.025)
    tmp.pred$uprq65 <- apply(q65, 1, quantile, probs = 0.975)
    tmp.pred$lwr25q65 <- apply(q65, 1, quantile, probs = 0.25)
    tmp.pred$upr75q65 <- apply(q65, 1, quantile, probs = 0.75)
    
    ## get size of CI for plotting
    tmp.pred$diffCI <- tmp.pred$upr - tmp.pred$lwr
    tmp.pred$diffCIq50 <- tmp.pred$uprq50 - tmp.pred$lwrq50
    tmp.pred$diffCIq65 <- tmp.pred$uprq65 - tmp.pred$lwrq65
    
    df.t <- cbind(df.t, tmp.pred[,colnames(tmp.pred)[!colnames(tmp.pred) %in% colnames(df.t)]]) ## check which rows to delete
    df.t$diff <- df.t$PE - df.t$pred_mod
    df.t$diffq50 <- df.t$PE - df.t$q50
    df.t$diffq65 <- df.t$PE - df.t$q65
    
    write.csv(df.t, paste0("./output_model/data_test_m5_", s, "_", m, ".csv"), row.names = F)
    
    print(ggplot(df.t) + 
      geom_abline(slope = 1, lty = "dashed", color = "blue") +
      geom_point(aes(x = PE, y = q50), alpha = 0.2, pch = 1) +
      geom_linerange(aes(x = PE, ymin = lwrq50, ymax = uprq50), lwd = 0.5, alpha = 0.1))
    
    print(ggplot(df.t) + 
      stat_halfeye(aes(x = Individual, y = diff)))
    print(ggplot(df.t) + 
      stat_halfeye(aes(x = Individual, y = diffq50)))
    print(ggplot(df.t) + 
      stat_halfeye(aes(x = Individual, y = diffq65)))
    
    print(ggplot(df.t) + 
      stat_halfeye(aes(x = Sc, y = diff)))
    print(ggplot(df.t) + 
      stat_halfeye(aes(x = Sc, y = diffq50)))
    print(ggplot(df.t) + 
      stat_halfeye(aes(x = Sc, y = diffq65)))
    
    dev.off()
  
  }
  write.csv(df.pred, paste0("./output_model/model-predictions_m5_", m, ".csv"), row.names = F)
  write.csv(df.pred2, paste0("./output_model/model-predictions2_m5_", m, ".csv"), row.names = F)
  write.csv(df.sim, paste0("./output_model/model-simulations_m5_", m, ".csv"), row.names = F)
  write.csv(df.coef, paste0("./output_model/model-coefficients_m5_", m, ".csv"), row.names = F)
}



#############
#############
cond <- predict(model, newdata = test_data, type = "response", re.form = NULL, cov.fit = T, se.fit = T)
cond2 <- predict(model, newdata = test_data, type = "link", re.form = NULL)
disp <- predict(model, newdata = test_data, type = "disp")
hist(cond)
hist(cond2)
hist(disp)
plot(x = cond, y = disp)
hist(rlnorm(1000, mean(cond), mean(disp)))
hist(rlnorm(1000, mean(log(cond)), mean(log(disp))))

m <- cond
v <- disp^2
# phi <- sqrt(v+m^2)
# mu <- log(m^2/phi)
sigma <- sqrt(log(1+v/m^2))
mu <- log(m)-0.5*sigma^2 # identisch zu oberem mu

test_data$sigma <- sigma
test_data$mu <- mu

set.seed(3171)
test_data <- test_data %>% rowwise() %>%
  mutate(pred1 = rlnorm(1, mu, sigma),
         pred2 = rlnorm(1, mu, sigma),
         pred3 = rlnorm(1, mu, sigma))
hist(test_data$pred1)
hist(test_data$pred2)
hist(test_data$pred3)
hist(test_data$PE)

ggplot(test_data) + geom_point(aes(x = PE, y = PE-pred1), color = "darkred", alpha = 0.5, size = 2, pch = 1) + 
  geom_point(aes(x = PE, y = PE-pred2), color = "lightblue", alpha = 0.5, size = 2, pch = 1) + 
  geom_point(aes(x = PE, y = PE-pred3), color = "grey", alpha = 0.5, size = 2, pch = 1)

bsim <- matrix(nrow = nrow(test_data), ncol = 1000)
for(i in 1:nrow(test_data)) {
  bsim[i, ] <- rlnorm(1000, test_data$mu[i], test_data$sigma[i])
}
hist(rowMeans(bsim))
hist(bsim[,1])
hist(bsim[,2])
hist(bsim[,3])

s <- "maisD"; m <- "direct.in"
mod <- readRDS(paste0("./output_model/model_m5_", s, "_", m, ".RDS"))
1 - sum((mod$frame$PE - predict(mod, type = "response"))^2)/
  sum((mod$frame$PE - mean(mod$frame$PE))^2)
