#-##############################################-#
# 1. Source the support file for this routine ----
#-##############################################-#
rm(list = ls()) # OPTIONAL: Execute to remove previous objects from the global environment.

# Activate the support 'functions file' for this routine:
source("./R_model/Linear modelling workflow_support functions.R") 
library(sf)
library(ggeffects) # to predict gams (in fast way)
library(parallel)
options(mc.cores = detectCores())


# Import data frame (example only. Any other import function of R can be used).
dsn <- paste0("../data/data_cali/savedFiles/Data_cali_density.gpkg")
dat <- st_read(layer = 'density', dsn = dsn)
dat <-  dat %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
dat.m <- st_read(layer = 'density_mean', dsn = dsn)
dat.m <-  dat.m %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

## change col names and merge
dat$pos <- "raw"
dat.m$pos <- "mean"

dat.m$PE <- dat.m$PE.m

dat <- rbind(dat[, c("site", "X_time", "Station.Count", "Antenna.Count", "Weight", 
                   "Signal.max", "Individual", "detR", "ant", "meth", "pos", "AcSc",
                   "lon", "lat", "lon.true", "lat.true", "PE", "dens")],
            dat.m[, c("site", "X_time", "Station.Count", "Antenna.Count", "Weight", 
                   "Signal.max", "Individual", "detR", "ant", "meth", "pos", "AcSc",
                   "lon", "lat", "lon.true", "lat.true", "PE", "dens")])
rm(dat.m)

# Carefully check data structure, column names and vector classes. Change them as needed.
str(dat)

dat$detR <- gsub('m', '', dat$detR)
dat$detR <- as.numeric(dat$detR)

## PE as log (then it's normal)
dat$PElog <- log10(dat$PE)

## add tag height
dat$tagheight_m[dat$Individual %in% c("TT090C", "TT241D")] <- 1.5
dat$tagheight_m[dat$Individual %in% c("TT163C", "TT298D", "TT164D", "TT014D")] <- 1.0
dat$tagheight_m[dat$Individual %in% c("TT240C", "TT090D")] <- 0.5
dat$tagheight_m[dat$Individual %in% c("TT298C")] <- 2.0


## create datasets for model 1-3 and maisC and D
## model 1 (ab_ql to get best detR)
df1C <- dat[dat$site == "maisC" & dat$meth == "ab_ql" & dat$pos == "raw" & dat$AcSc == "ac1-sc1",]
df1D <- dat[dat$site == "maisD" & dat$meth == "ab_ql" & dat$pos == "raw" & dat$AcSc == "ac1-sc1",]

## model 2 (get Ac Sc threshold)
df2C <- df1C[df1C$detR == 900, ]
df2D <- df1D[df1D$detR == 800, ]
df22C <- dat[dat$site == "maisC" & dat$meth == "in_ql" & dat$pos == "raw" & dat$AcSc == "ac4-sc2",]
df22D <- dat[dat$site == "maisD" & dat$meth == "in_ql" & dat$pos == "raw" & dat$AcSc == "ac4-sc2",]
df23C <- dat[dat$site == "maisC" & dat$meth == "ab_ml" & dat$pos == "raw" & dat$AcSc == "ac1-sc1",]
df24C <- dat[dat$site == "maisC" & dat$meth == "ml_ml" & dat$pos == "raw" & dat$AcSc == "ac2-sc2",]

## model 3 (mean or raw position)
df3C <- dat[dat$site == "maisC" & dat$meth == "ab_ql" & 
              dat$AcSc == "ac3-sc1" & 
              dat$detR == 900 ,]
df3D <- dat[dat$site == "maisD" & dat$meth == "ab_ql" & 
              dat$AcSc == "ac3-sc1" & 
              dat$detR == 800 ,]
df32C <- dat[dat$site == "maisC" & dat$meth == "in_ql" & dat$AcSc == "ac4-sc2",]
df32D <- dat[dat$site == "maisD" & dat$meth == "in_ql" & dat$AcSc == "ac4-sc2",]
df33C <- dat[dat$site == "maisC" & dat$meth == "ab_ml" & dat$AcSc == "ac2-sc2",]
df34C <- dat[dat$site == "maisC" & dat$meth == "ml_ml" & dat$AcSc == "ac2-sc2",]

## model 4 (best method)
df4C <- do.call(rbind, list(df3C, df32C, df33C, df34C))
df4C <- df4C[df4C$pos == "mean",]
df4D <- do.call(rbind, list(df3D, df32D))
df4D <- df4D[df4D$pos == "mean",]

## model 5 (position accuracy for best method)
df5C <- df4C[df4C$meth == "ab_ql",]
df5D <- df4D[df4D$meth == "ab_ql",]

## define dataset
df <- df5C

## get factors
df$Individual <- as.factor(df$Individual)
df$ant <- as.factor(df$ant)
df$meth <- as.factor(df$meth)
df$pos <- as.factor(df$pos)
df$AcSc <- as.factor(df$AcSc)

df$detR <- factor(df$detR, ordered = T)
df$Sc <- factor(df$Station.Count, ordered = T)
df$Ac <- factor(df$Antenna.Count, ordered = T)
df$tagheight_m <- as.factor(df$tagheight_m)

## Model 1 ---------------------------------------------------------------------
# get best detR

var_resp <- "PE"   

var_fac <- "detR"
var_num <- NULL

var_fac <- c("Ac", "Sc")
var_num <- c("Signal.max")

var_fac <- NULL
var_num <- c("Signal.max", "Station.Count")

var_fac <- c("Sc")
var_num <- c("Signal.max")

var_fac <- c("pos")
var_num <- NULL

var_fac <- c("meth")
var_num <- NULL

var_fac <- NULL
var_num <- c("Antenna.Count", "Station.Count", "Weight", "Signal.max", "dens")

var_rand <- "tagheight_m"

df <- remove_NAs(data = df, variables = c(var_num, var_fac, var_rand, var_resp))

## normal sample
samp <- 20000
df <- df[sample(c(rep(TRUE, samp), rep(FALSE, nrow(df)-samp)), nrow(df) ,replace = F),]

## sample timestamps present in all methods (for mod 4)
df <- df %>% group_by(X_time) %>%
  mutate(Nmeth = length(unique(meth))) %>%
  ungroup()

df <- df[df$Nmeth == nlevels(df$meth),]

df <- cbind(df, z_transform(data = df, predictors = var_num))
corvif(data = df, variables = c(var_num, var_fac))

# mod1 <- glmmTMB(PE ~ detR + (1|tagheight_m),
#                data = df,                          # Name of dataframe.
#                family = Gamma(link = "log"))     # Template: Replace with your distribution family and link.

# mod2 <- glmmTMB(PE ~ 
#                   # Sc + Ac + 
#                   # Station.Count_z*Antenna.Count_z + 
#                   # Station.Count_Z +
#                   Sc + 
#                   Signal.max_z + (1|tagheight_m),
#                 data = df,                          # Name of dataframe.
#                 family = Gamma(link = "log"))     # Template: Replace with your distribution family and link.

# mod3 <- glmmTMB(PE ~ 
#                   pos + (1|tagheight_m),
#                 data = df,                          # Name of dataframe.
#                 family = Gamma(link = "log"))     # Template: Replace with your distribution family and link.

# mod4 <- glmmTMB(PE ~ 
#                   meth + (1|tagheight_m),
#                 data = df,                          # Name of dataframe.
#                 family = Gamma(link = "log"))     # Template: Replace with your distribution family and link.

mod5 <- glmmTMB(PE ~ 
                  Station.Count_z*Antenna.Count_z*dens_z +
                  Signal.max_z + Weight_z +
                  (1|tagheight_m),
                data = df,                          # Name of dataframe.
                family = Gamma(link = "log"))     # full model has convergence issues (all interactions)

# library(brms)
# library(cmdstanr)
# 
# mod5 <- brm(PE ~ 
#                   Station.Count_z*Antenna.Count_z +
#                   dens_z +
#                   Signal.max_z + Weight_z +
#                   (1|tagheight_m),
#                 data = df,                          # Name of dataframe.
#             backend = "cmdstanr", 
#             threads = threading(3),
#             cores = 4,
#             
#                 family = Gamma(link = "log"))


mod <- mod5

# This function plots the overall distribution of model residuals.
residual_plots(data = df,
               modelTMB = mod,
               response = "PE") # Name of response variable.

# This function plots residuals against all possible predictors specified in 1.
residual_plots_predictors(data = df,
                          modelTMB = mod,
                          predictors = c(var_num, var_fac)) # Name of all fixed predictors in the dataframe.

# Graphical displays of all 2-way predictor-response relationships.
if(length(c(var_num, var_fac)) >= 2) {
  residual_plots_interactions(data = df, 
                            modelTMB = mod, 
                            predictors = c(var_num, var_fac))      # All fixed predictors in the dataframe.

}

# These plots can reveal if any random intercept term additionally needs random slopes over a specific fixed predictor.
residual_plots_random_slopes(data = df,                       # Name of dataframe.
                             modelTMB = mod,                  # Name of model.
                             fixed_eff = c(var_num, var_fac), # Name of fixed factors meeting criterion described above.
                             random_eff = var_rand)           # Name of random factors meeting criterion described above.


# Compare the variance of the observed data with the variance distribution in model-simulated dataframes:
dispersion_simulation(data = df, 
                      modelTMB = mod,
                      response = var_resp,
                      predictor = var_fac[1], # (optional, but recommended) - A SINGLE factor predictor to split this plot. Use a different index number as needed.
                      n.sim = 500)            # Number of simulations. Set to >2000 for final model assessment.

mean_simulation(data = df, 
                      modelTMB = mod,
                      response = var_resp,
                      predictor = var_fac[1], # (optional, but recommended) - A SINGLE factor predictor to split this plot. Use a different index number as needed.
                      n.sim = 500)            # Number of simulations. Set to >2000 for final model assessment.

sd_simulation(data = df, 
                modelTMB = mod,
                response = var_resp,
                #predictor = var_fac[1], # (optional, but recommended) - A SINGLE factor predictor to split this plot. Use a different index number as needed.
                n.sim = 500)            # Number of simulations. Set to >2000 for final model assessment.

# The function compares the observed raw data distribution with that observed in model-simulated dataframes.
ppcheck_fun(data = df,
            modelTMB = mod,
            response = var_resp,
            #predictor = var_fac[1],   # (optional, but recommended) - A SINGLE factor predictor to split this plot.
            n.sim = 500)                # Number of simulations, set to >2000 for final model assessment.


# Extract coefficient estimates and 95% compatibility intervals:
comp_int(modelTMB = mod,
           ci_range = 0.95,       # Compatibility interval range.
           effects = "all",       # Returns parameters for fixed effects ("fixed"), random effects ("random"), or both ("all").
           component = "all")     # Returns parameters for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").


# Extract marginal and conditional R-squared value.
r2(model = mod) # Name of model
# NOTE: this works for (G)LMMs, only, i.e., for models that include a random component.

if(!is.null(var_num) & !is.null(var_fac)) {
  for(i in var_num) {
    
    print(slope_estimates(data = df,
                          modelTMB = mod,
                          num_predictor = i,  # Name of one numeric predictor for which to estimate slopes. Change index as needed.
                          fac_predictors = var_fac[1]) # # Name of one (or two!) factor predictor(s). Slope estimates will be given per level of this factors. Change index as needed. 
          
    )
  }
}

#==========================================================#
# * 7.3 Estimation of differences between factor levels ---- 
#==========================================================#
# Extract pairwise comparisons among levels of specified factor predictors (means and 95% compatibility intervals).
if(!is.null(var_fac)) {
  pairwise_comparisons(data = df,
                       modelTMB = mod,
                       predictors = var_fac[1], # Name of one (or more!) factor predictor(s).
                       component = "cond",       # Returns differences for the conditional model ("cond"), or the zero-inflation part of the model ("zi").
                       dispScale = "response",  # "response" returns absolute differences for link identity models, and ratio or odds.ratios for log or logit link models 
                       # "link" returns estimates on the link scale
                       contrasts = "all")       # "all" returns all pairwise comparisons, "within" only the comparisons among the levels of the first factor specified within each level of the second one
  
}

newdat <- expand.grid("detR" = factor(unique(df$detR), ordered = T),
                      "pos" = levels(df$pos),
                      "meth" = levels(df$meth),
                      # "AcSc" = levels(df$AcSc),
                      "dens_z" = c(-1, 0, 1),
                      "Station.Count" = unique(df$Station.Count),
                      "Antenna.Count" = unique(df$Antenna.Count),
                      # "Sc" = factor(unique(df$Sc), ordered = T),
                      # "Ac" = factor(unique(df$Ac), ordered = T),
                      "Signal.max_z" = mean(df$Signal.max_z),
                      "Weight_z" = mean(df$Weight_z),
                      "tagheight_m" = as.factor(unique(df$tagheight_m)))

newdat$Station.Count_z <- (newdat$Station.Count-mean(df$Station.Count))/sd(df$Station.Count)
newdat$Antenna.Count_z <- (newdat$Antenna.Count-mean(df$Antenna.Count))/sd(df$Antenna.Count)

## reduce to plausible values for Sc and Ac
# newdat <- newdat[as.numeric(newdat$Sc) <= as.numeric(newdat$Ac) &
#                    as.numeric(newdat$Sc) >= as.numeric(newdat$Ac)/4,]
newdat <- newdat[newdat$Station.Count <= newdat$Antenna.Count & 
                   newdat$Station.Count >= newdat$Antenna.Count/4,]

for(i in 1:nrow(newdat)) {
  Ac <- newdat$Antenna.Count[i]
  Sc <- newdat$Station.Count[i]
  newdat$dens_z[i] <- mean(df$dens_z[df$Antenna.Count == Ac & df$Station.Count == Sc])
  newdat$Signal.max_z[i] <- mean(df$Signal.max_z[df$Antenna.Count == Ac & df$Station.Count == Sc])
  newdat$Weight_z[i] <- mean(df$Weight_z[df$Antenna.Count == Ac & df$Station.Count == Sc])
  
}

newdat <- newdat[!is.nan(newdat$dens_z),] # remove NaNs due to non-present combinations

# pred <- predict(mod, newdat, type = "response", re.form = NA, se.fit = T)
# newdat$median <- pred[["fit"]]
# newdat$se <- pred[["se.fit"]]

mod_predictionsN <- post_predictN(data = df,
                                 modelTMB = mod,
                                 # newdat = dplyr::select(newdat, -c(Station.Count)),
                                 newdat = dplyr::select(newdat, -c(Station.Count, Antenna.Count)),
                                 # newdat  = newdat, # Name of predictors specified in 8.1.
                                 # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                                 component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").
ggplot(mod_predictionsN) + 
  geom_pointrange(aes(x = meth, ymin = lower, ymax = upper, y = median)) +
  ylab("position error") +
  theme_light() 

ggplot(mod_predictionsN) + 
  geom_pointrange(aes(x = pos, ymin = lower, ymax = upper, y = median)) +
  ylab("position error") +
  theme_light() 

ggplot(mod_predictionsN) + 
  geom_pointrange(aes(x = detR, ymin = lower, ymax = upper, y = median)) +
  ylab("position error") +
  theme_light() 

ggplot(mod_predictionsN) + 
  geom_pointrange(aes(x = Sc, ymin = lower, ymax = upper, y = median)) +
  ylab("position error") +
  theme_light() 

ggplot(mod_predictionsN) + 
 # geom_point(aes(x = Sc, y = Ac), data = df, pch = 1, size = 3) +
  # geom_point(aes(x = Sc, y = Ac, color = median, size = SD)) +
  geom_point(aes(x = Station.Count, y = Antenna.Count, color = median, size = SD)) +
  scale_color_viridis_c(option = "rocket", limits = c(NA, NA)) +
  facet_wrap(~as.factor(dens)) +
  theme_light()  


ggplot(newdat) + 
  # geom_point(aes(x = Sc, y = Ac), data = df, pch = 1, size = 3) +
  # geom_point(aes(x = Sc, y = Ac, color = median, size = SD)) +
  # geom_point(aes(x = Station.Count, y = Antenna.Count, color = median, size = se)) +
  geom_pointrange(aes(x = Station.Count, y = median, ymin = median - 1.96*se, ymax = median + 1.96*se)) +
  scale_color_viridis_c(option = "rocket") +
  theme_light() 

df.sum <- df %>% group_by(Ac, Sc) %>%
  summarize(mPE = mean(PE), 
            nP = n(), .groups = "drop")

ggplot(df.sum) + 
  geom_point(aes(x = Sc, y = Ac, color = mPE, size = nP)) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  theme_light() 


pairwise_comparisons(data = df,
                     modelTMB = mod,
                     predictors = var_fac[1], # Name of one (or more!) factor predictor(s).
                     component = "cond",       # Returns differences for the conditional model ("cond"), or the zero-inflation part of the model ("zi").
                     dispScale = "response",  # "response" returns absolute differences for link identity models, and ratio or odds.ratios for log or logit link models 
                     # "link" returns estimates on the link scale
                     contrasts = "all")       # "all" returns all pairwise comparisons, "within" only the comparisons among the levels of the first factor specified within each level of the second one

ggplot(df) + geom_point(aes(x = Weight, y = PE), alpha = 0.5, pch = 1) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  facet_wrap(~Station.Count) +
  theme_light()

ggplot(df) + geom_point(aes(x = Weight, y = PE), alpha = 0.5, pch = 1) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  facet_wrap(~Antenna.Count) +
  theme_light()

ggplot(df) + geom_point(aes(x = Signal.max, y = PE), alpha = 0.5, pch = 1) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  facet_wrap(~Station.Count) +
  theme_light()

ggplot(df) + geom_point(aes(x = Signal.max, y = PE), alpha = 0.5, pch = 1) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  facet_wrap(~Antenna.Count) +
  theme_light()

ggplot(df) + geom_point(aes(x = dens, y = PE), alpha = 0.5, pch = 1) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  facet_wrap(~Station.Count) +
  theme_light()

ggplot(df) + geom_point(aes(x = dens, y = PE), alpha = 0.5, pch = 1) +
  scale_color_viridis_c(option = "rocket", limits = c(0, 150), na.value = "#FAEBDDFF") +
  facet_wrap(~Antenna.Count) +
  theme_light()

## bootstrapping data instead of mean predictions
# Bootstrap predictions (using simulate repeatedly)
bootstrap_responses <- replicate(1000, unlist(simulate(mod, newdata = newdat)))

# Compute prediction intervals
tmp <- apply(bootstrap_responses, 1, quantile, probs = c(0.025, 0.5, 0.975))
