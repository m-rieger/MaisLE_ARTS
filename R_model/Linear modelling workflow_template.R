#-##############################################-#
# 1. Source the support file for this routine ----
#-##############################################-#
rm(list = ls()) # OPTIONAL: Execute to remove previous objects from the global environment.

# Activate the support 'functions file' for this routine:
source("./R_model/Linear modelling workflow_support functions.R") 
library(sf)
library(ggeffects) # to predict gams (in fast way)

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
df1C <- dat[dat$site == "maisC",]
df1C <- df1C[df1C$meth == "ab_ql",]

df1D <- dat[dat$site == "maisD",]
df1D <- df1D[df1D$meth == "ab_ql",]

## model 2 (best ab_ql and other meth)
df2C <- dat[dat$site == "maisC",]
df2C <- df2C[df2C$detR %in% c(NA, 800),]

df2D <- dat[dat$site == "maisD",]
df2D <- df2D[df2D$detR %in% c(NA, 800),]

## model 3 (optimal pos. error estimation)
df3C <- dat[dat$site == "maisC",]
df3C <- df3C[df3C$detR == 800,]
df3C <- df3C[df3C$meth == "ab_ql",]
df3C <- df3C[df3C$pos == "mean",]

df3D <- dat[dat$site == "maisD",]
df3D <- df3D[df3D$detR == 800,]
df3D <- df3D[df3D$meth == "ab_ql",]
df3D <- df3D[df3D$pos == "mean",]

## define dataset
df <- df1D

## get factors
df$Individual <- as.factor(df$Individual)
df$ant <- as.factor(df$ant)
df$meth <- as.factor(df$meth)
df$pos <- as.factor(df$pos)
df$AcSc <- as.factor(df$AcSc)

df$detR <- as.factor(df$detR)

#-##############################-#
# 3. Definition of variables  ----
#-##############################-#
# In step 3.1 through 3.4, the names of all variables (= data frame columns) that are potentially relevant 
# for your analysis are stored as STRINGS for further use:
#================================#
# * 3.1 ONE RESPONSE variable ----
#================================#
var_resp <- "PE"                     # Replace with the name of the response variable.


#=============================================================================#
# * 3.2 Fixed predictors: quantitative and categorical predictor variables ----
#=============================================================================#
# FACTOR PREDICTOR variable(s)
var_fac <- c("meth", "pos", "AcSc", "detR")                              # assign default NA if missing.
var_fac2 <- c("pos", "AcSc", "detR")                              # assign default NA if missing.
# var_fac <- c("fac_1", "fac_2", ... )   # Template: replace entries with the name(s) of your factor predictor(s).
                                           # Assure all these are factors with at least two levels.

# NUMERIC or INTEGER PREDICTOR variable(s) 
var_num <- c("dens", "Station.Count", "Antenna.Count", "Signal.max", "tagheight_m", "Weight")                              # assign default NA if missing.
var_num2 <- c("dens", "Station.Count", "Antenna.Count", "Signal.max", "Weight")                              # assign default NA if missing.
# var_num <- c("num_1", "num_2", ... )   # Template: replace entries with the name(s) of your numeric predictor(s).
                                           # Assure all these are numeric or integer.

#==============================================================#
# * 3.3 Random predictors: dependency structure of the data ----
#==============================================================#
# RANDOM term(s)
var_rand <- c("Individual") # Individual                        # assign default NA if missing.
  # var_rand <- c("rand_1", "rand_2", ...) # Template: replace entries with the name(s) of your factor random term(s).
                                           # Assure all these are factors with at least five levels.

#==============================================#
# * 3.4 Temporal and spatial data structure ----
#==============================================#
# The variables specified here will be used to check for temporal and/or spatial autocorrelation in model residuals (6.3).

# NUMERIC or INTEGER TIME variable that specifies temporal structure.
  # First, enter ONE variable that has information on temporal sequence/time:
var_time <- NA                          # assign default NA if missing.
    # var_time <- "daytime"               # Template: replace entry with the name of your temporal variable.
  
  # Second, if these temporal data are structured, add the grouping variable below.
  # (example: multiple parallel time-series that are split by experimental blocks)
  var_time_groups <- NA                   # assign NA if missing.
    # var_time_groups <- "daytime.grouping" # Template: replace entry with the name of your temporal grouping variable.

# NUMERIC or INTEGER COORDINATES (x, y) that specify spatial structure.
  var_space <- NA                         # assign NA if missing.
    # var_space <- c("x_coord","y_coord") # Template: replace entry with the name of your variables that contain spatial coordinates.

  
#=========================#
# * 3.5 Missing values ----
#=========================#
# We prune the dataset to complete cases for all variables that are considered for analysis.
df.pr <- remove_NAs(data = df, variables = c(var_num, var_fac, var_rand, var_resp))
  # NOTE-1: Data pruning should be kept to the minimum needed to allow models to run. 
  #         Therefore, rerun 1. after a final model formulation has been identified. This may allow to rescue some observations.
  # NOTE-2: A key assumption is that missing data (and the reasons for data to be missing) are randomly distributed across the dataset. 
  #         Confirm by inspecting the removed data rows after creating the dataset df.NAs 
            df.NAs <- anti_join(df,  df.pr) # This object contains all removed data rows.
  
# Keep only complete cases in the dataset for further analysis: 
df <- df.pr
table(df$pos, df$AcSc, df$Individual)

## for TESTING: only part of data
samp <- 20000
df <- df[sample(c(rep(TRUE, samp), rep(FALSE, nrow(df)-samp)), nrow(df) ,replace = F),]

#-###############################-#
# 4. Raw data exploration ----
#-###############################-#
# This step helps to identify patterns that may require consideration for model formulation.
# NOTE: Further checks can be found in the residual plot analysis (section 6.1), in particular for missing predictors or interactions.

#=========================#
# * 4.1 Extreme values ----
#=========================#
# We graphically inspect variables for extreme values.

#---------------------------------------------#
# ** 4.1.1 Extremes in NUMERICAL variables ---- 
#---------------------------------------------#
dotplot_num(data = df, variables = c(var_num, var_resp)) # Dotplots for all numerical predictors plus the response variable.

# What should I look for?
# >> Is any observation CLEARLY separated from the core distribution? Such observersation may represent 'implausible extremes'

# Resolving implausible extreme values: 
# >> Check where such values originate from. If you can trace them to objective (!) typing errors, correct these values.
# >> If extremes are in the response variable: choose an adequate distribution family.
# >> If extremes are in predictor variables: 
#    1. Use model assessment to check if these predictor values cause concern.
#    2. If so, consider data transformation to mitigate the issue.
# >> If extremes cannot be modeled appropriately, consider reporting effect estimates and their SE with and without these extremes.


#------------------------------------------#  
# ** 4.1.2 Extremes in FACTOR variables ----      
#------------------------------------------#
barplot_fac(data = df, variables = c(var_fac, var_rand)) # Barplots for all factor predictors and random terms.

# What should I look for?
# >> Are observations roughly balanced between grouping levels (conforming to a "balanced design")?

# Resolving extreme unbalance: 
# >> If some factor levels arise from typing errors: Correct them.
# >> If there is extreme imbalance in sample size among factor levels, 
#    accept that effect estimates for those with poor replication will associate with large uncertainty.
#    Alternatively, consider pooling levels with very few replicates, but only if this remains biologically meaningful. 

  
#================================#
# * 4.2 Predictor collinearity ----
#================================#
# We graphically and numerically inspect variables for predictor collinearity.

#-------------------------------------------------------------#
# ** 4.2.1: Graphical inspection for predictor collinearity ----
#-------------------------------------------------------------#
# We graphically inspect predictor collinearity for all pairwise combinations of numeric and/or factor predictors.

# *** 4.2.1.1: NUMERIC predictors: Scatterplots ----
# Skip this step if you have < 2 numeric predictors
coll_num(data = df, predictors = var_num) # Pairwise scatterplots for all numeric predictors.
# What should I look for?
# >> Data should distribute rather homogeneously.


# *** 4.2.1.2: FACTOR against NUMERIC predictors: Swarm plots ----
# Skip this step if your predictors are either all numeric, OR all factor.
coll_num_fac(data = df, predictors = c(var_num, var_fac))  # Swarm boxplots for all numeric against all factor predictors.
# What should I look for?
# >> Numeric values on the y axis should distribute roughly homogeneously across levels on the x-axis.
 

# *** 4.2.1.3: FACTOR predictors: Mosaic plots ----   
# Skip this step if you have < 2 factor predictors.
coll_fac(data = df, predictors = var_fac)  # Pairwise mosaic plots for all factor predictors.
# What should I look for?
# >> All combinations of predictor levels should have roughly similar sample sizes.


#------------------------------------------------#
# ** 4.2.2: Variance Inflation Factors (VIFs) ----
#------------------------------------------------#
corvif(data = df, variables = c(var_num2, var_fac2))  # VIF are calculated for all fixed predictors potentially included in the model.
# Derived from Zuur, Ieno & Elphick (2010).

# What should I look for?
# >> Check column "GVIF": Entirely independent predictors yield GVIF = 1, while larger GVIF values indicate 
#   increasing predictor collinearity (see article section 4.2.2 for details)

#==========================================#
#* 4.3 Predictor-response relationships ----
#==========================================#
# The function generates plots of the response variable against each of the potential predictor variables.
relat_single(data = df,
             response = var_resp,                        # Name of response variable
             # y.limits = c(min,max),                    # (optional) - Limits of the y-axis. Replace min and max with your limits
             # y.breaks = c(value_1,value_2,value_3,...),# (optional) - Breaks of the y-axis. Replace with your values
             predictors = c(var_num, var_fac)            # Name of fixed predictors
             )  
# NOTE: Up to 9 plots are shown in the plot window... 
#       If > 9 plots are present, they are saved in your working directory in a file 
#       named "Predictor_response_relationships.pdf"

# What should I look for?
# >> Identify relevant relationships with the response.


#===============================#
# * 4.4 Response distribution ----
#===============================#
# Visualise the distribution of your response variable.
distr_response(data = df, response = var_resp) 
  # What should I look for?
  # >> The observed distribution of the response confirms - at least roughly - your initial expectation.
  #    If not, adjust the distribution family for your initial model accordingly. 


#-#######################-#
# 5. Model formulation ----
#-#######################-#
# Based on the conclusions from data exploration, formulate your preliminary model.

#===============================================#
# * 5.1 Standardisation of numeric predictors ----
#===============================================#
# This step is optional, but standardisation is recommended.
# The code adds z-transformed numeric covariates to the dataframe. We recommend to use these for model formulation.
df <- cbind(df, z_transform(data = df,         # z-transformation ...
                        predictors = var_num)) # ... for all numeric predictor variables.

    # Users who have already standardised their predictors have two options:
    # 1. Skip this step. and accept that final results plots will display standardised predictor values.
    # 2. Start with predictors on the raw scale, and use our function for scaling to allow correct final plotting.


#===============================#
# * 5.2 Model implementation ----
#===============================#
library(parallel)
options(mc.cores = detectCores())

#---------------------------------------#  
# ** 5.2.1 Initial model formulation ----
#---------------------------------------#  
# Implement the initial model.
# See article section 5.2 for guidance to model formulation.


mod1 <- glmmTMB(PE ~
                 detR + detR:pos +
                 Weight_z + Weight_z:pos +
                 Antenna.Count_z + Antenna.Count_z:pos +
                 Station.Count_z + Station.Count_z:pos +
                 dens_z + dens_z:pos +
                 pos +
                 AcSc +
                 (1|Individual),

               data = df,                          # Name of dataframe.
               family = Gamma(link = "log"))     # Template: Replace with your distribution family and link.
mod <- mod1
# mod2 <- glmmTMB(PE ~ 
#                  meth + meth:pos +
#                  dens_z + dens_z:meth +
#                  Station.Count_z + Station.Count_z:meth +
#                  Signal.max_z + Signal.max_z:meth +
#                  # tagheight_m_z +
#                  pos + 
#                  (1|Individual),         
#                
#                data = df,                          # Name of dataframe.
#                family = Gamma(link = "log"))     # Template: Replace with your distribution family and link.


# mod.ln <- glmmTMB(PE ~ 
#                  meth + meth:pos +
#                  dens_z + dens_z:meth +
#                  Station.Count_z + Station.Count_z:meth +
#                  Signal.max_z + Signal.max_z:meth +
#                  # tagheight_m_z +
#                  pos + 
#                  (1|Individual),         
#                
#                data = df,                          # Name of dataframe.
#                family = lognormal(link = "log"))     # Template: Replace with your distribution family and link.


# This formulation may need refinement (section 5.2.2 in manuscript) after model assessment (6.).


#-######################-#
# 6. Model assessment ----  
#-######################-#

# Before inspecting model results, careful assessment is required.

#========================#
#* 6.1 Residual checks ----
#========================#
# This step helps with checking for 'residual patterns'.

#-----------------------------------------#
# ** 6.1.1: Distribution of residuals -----
#-----------------------------------------#
# This function plots the overall distribution of model residuals.
residual_plots(data = df,
               modelTMB = mod,
               response = var_resp) # Name of response variable.

# What should I look for? 
# >> QQ-plot for residuals: Points should closely follow the diagonal reference line.
# >> Residuals against fitted: Points should homogeneously scatter without any pattern.
# >> QQ-plot for random intercepts, ONLY produced for (G)LMMs: Points should closely follow the diagonal reference line.

# Resolving violations:
# >> Seek a distribution family and/or a link function that better capture the observed residual distribution.
# >> The model may lack an important covariate (-> 6.1.2), or an informative interaction term (-> 6.1.3).
# >> The model may lack relevant non-linear (polynomial) terms (-> 6.1.2).
# >> The model could be over- or underdispersed (-> 6.2.1).
# >> The model may suffer from zero-inflation (-> 6.2.2).


#---------------------------------------------------#
# ** 6.1.2: Residuals against possible PREDICTORS ----
#---------------------------------------------------#
# This function plots residuals against all possible predictors specified in 1.
residual_plots_predictors(data = df,
                          modelTMB = mod,
                          predictors = c(var_num, var_fac)) # Name of all fixed predictors in the dataframe.
# NOTE: For this check, we recommend that var_num and var_fac ALSO contain available variable(s) 
# that are currently NOT part of the model. 

# What should I look for? 
# >> Points should homogeneously scatter without obvious pattern.

# Resolving violations:
# >> Seek a distribution family and/or link function that better capture the observed residual distribution.
# >> The model may lack an important covariate (shown here: covariates that show residual patterns).
#    or an informative interaction term (-> 6.1.3).
# >> The model may lack relevant non-linear (polynomial) terms (shown here: curvilinear patterns in residuals across a given covariate).
# >> The model could be over- or underdispersed (-> 6.2.1).
# >> The model may suffer from zero-inflation (-> 6.2.2).
# >> Add a dispersion formula to the model to explicitly integrate heterogeneous variance across predictor values / levels.


#------------------------------------------------------------------------------#  
# ** 6.1.3: Residuals against possible TWO-WAY INTERACTIONS AMONG PREDICTORS ----
#------------------------------------------------------------------------------#  
# Graphical displays of all 2-way predictor-response relationships.
residual_plots_interactions(data = df, 
                            modelTMB = mod, 
                            predictors = c(var_num2, var_fac2))      # All fixed predictors in the dataframe.
# NOTE-1: In numeric vs. numeric display panels, smoothing lines are added as a visual aid to detect non-linear relationships.
# NOTE-2: Up to 9 plots are shown in the plot window. 
#         If > 9 plots are present, they are saved in your working directory as "Two-ways_interactions.pdf".

# What should I look for?   
# >> Points should homogeneously scatter without obvious pattern. 

# Resolving violations:
# >> Consider adding the relevant interaction term to the model if meaningful. 


#------------------------------------------------------------------------#
# ** 6.1.4: Residuals against PREDICTORS split by RANDOM FACTOR LEVELS ----
#------------------------------------------------------------------------#
# These plots can reveal if any random intercept term additionally needs random slopes over a specific fixed predictor.
  # NOTE: ONLY useful when MULTIPLE observations per combination of random and factor levels are present.
  #       Please only select those var_fac, var_num and/or var_rand that meet this criterion.
residual_plots_random_slopes(data = df,                       # Name of dataframe.
                             modelTMB = mod,                  # Name of model.
                             fixed_eff = c(var_num, var_fac), # Name of fixed factors meeting criterion described above.
                             random_eff = var_rand)           # Name of random factors meeting criterion described above.

# What should I look for?
# >> Do plot panels indicate diverging slopes (for numeric predictors),
#    or inconsistent effect directions (for factor predictors)?

# What if this is the case? 
# >> Add a random slope for that specific predictor to the model (see article section 5.2.1)


#=====================================#
# * 6.2 Posterior predictive checks ----
#=====================================# 
# Assess model performance with posterior predictive checks.
# These compare the observed data with the distribution of many, ideally >2000, 
# replicate raw dataframes simulated from the model. For initial checks a lower number of n.sim will do, though.

#-----------------------#
# ** 6.2.1 Dispersion ----
#-----------------------#
# Compare the variance of the observed data with the variance distribution in model-simulated dataframes:
dispersion_simulation(data = df, 
                      modelTMB = mod,
                      response = var_resp,
                      predictor = var_fac2[1], # (optional, but recommended) - A SINGLE factor predictor to split this plot. Use a different index number as needed.
                      n.sim = 500)            # Number of simulations. Set to >2000 for final model assessment.

# What should I look for?
# >> Observed variance should be roughly central within the simulated distributions.
# >> Check the function's feedback messages in the R console window.

# Resolving violations:
# >> Check if dispersion issues are connected to zero inflation (-> 6.2.2).
# >> Look for a distribution family and link that better capture the observed data dispersion.

mean_simulation(data = df, 
                      modelTMB = mod,
                      response = var_resp,
                      predictor = var_fac2[1], # (optional, but recommended) - A SINGLE factor predictor to split this plot. Use a different index number as needed.
                      n.sim = 500)            # Number of simulations. Set to >2000 for final model assessment.

sd_simulation(data = df, 
                modelTMB = mod,
                response = var_resp,
                predictor = var_fac2[1], # (optional, but recommended) - A SINGLE factor predictor to split this plot. Use a different index number as needed.
                n.sim = 500)            # Number of simulations. Set to >2000 for final model assessment.

#-------------------------------#
# ** 6.2.3 Data distribution ----
#-------------------------------# 
# The function compares the observed raw data distribution with that observed in model-simulated dataframes.
ppcheck_fun(data = df,
            modelTMB = mod,
            response = var_resp,
            predictor = var_fac2[1],   # (optional, but recommended) - A SINGLE factor predictor to split this plot.
            n.sim = 500)                # Number of simulations, set to >2000 for final model assessment.

# What should I look for?
# >> Ideally, the observed data distribution should be roughly central within the pattern 
#    of simulated data from the model.

# Resolving violations:
# >> Carefully re-evaluate the preceding steps of model assessment, and reformulate your model.

#-###########################################-#
# 7. Model results and parameter estimates ----
#-###########################################-#
# Perform this step only after the iteration of 3. and 4. results in a final model.

#=============================================================#
#* 7.1 Overall coefficient estimates and model performance ---- 
#=============================================================#
# Extract coefficient estimates and 95% compatibility intervals:
comp_int(modelTMB = mod,
           ci_range = 0.95,       # Compatibility interval range.
           effects = "all",       # Returns parameters for fixed effects ("fixed"), random effects ("random"), or both ("all").
           component = "all")     # Returns parameters for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").


# Extract marginal and conditional R-squared value.
r2(model = mod) # Name of model
# NOTE: this works for (G)LMMs, only, i.e., for models that include a random component.

mod.r <- update(mod, .~
                  detR_z + detR_z:pos +
                  Weight_z + Weight_z:pos +
                  Antenna.Count_z + Antenna.Count_z:pos +
                  Station.Count_z + Station.Count_z:pos +
                  dens_z + dens_z:pos +
                  pos +
                  (1|Individual))

mod.wi <- update(mod, .~
                  poly(detR_z, 2) + 
                  Weight_z +
                  Antenna.Count_z +
                  Station.Count_z +
                  dens_z + 
                  pos +
                  (1|Individual))

mod.f <- update(mod, .~
                  poly(detR_z, 2) + poly(detR_z, 2):pos +
                  poly(Weight_z, 2) + poly(Weight_z, 2):pos +
                  poly(Antenna.Count_z, 2) + poly(Antenna.Count_z, 2):pos +
                  poly(Station.Count_z, 2) + poly(Station.Count_z, 2):pos +
                  dens_z + dens_z:pos +
                  pos +
                  (1|Individual))

mod.s <- update(mod, .~
                  s(detR_z, k = 5) +
                  Weight_z + Weight_z:pos +
                  Antenna.Count_z + Antenna.Count_z:pos +
                  Station.Count_z + Station.Count_z:pos +
                  dens_z + dens_z:pos +
                  pos +
                  (1|Individual))
mod.0 <- update(mod, .~
                  poly(detR_z, 2) + poly(detR_z, 2):pos +
                  pos +
                  (1|Individual))
r2(mod); r2(mod.r); r2(mod.s); r2(mod.0); r2(mod.f); r2(mod.wi)

#=========================================#
#* 7.2 Estimation of regression slopes ---- 
#=========================================#
# Extract slope estimates, their 95% compatibility intervals, and an estimate of (proportional) change from simulated data.
# This function is useful to extract slope estimates within each level of a factor predictor that is also part of the model.

slope_estimates(data = df,
                modelTMB = mod,
                num_predictor = "dens_z",  # Name of one numeric predictor for which to estimate slopes. Change index as needed.
                fac_predictors = var_fac2[1]) # # Name of one (or two!) factor predictor(s). Slope estimates will be given per level of this factors. Change index as needed. 

slope_estimates(data = df,
                modelTMB = mod,
                num_predictor = "Antenna.Count_z",  # Name of one numeric predictor for which to estimate slopes. Change index as needed.
                fac_predictors = var_fac2[1]) # # Name of one (or two!) factor predictor(s). Slope estimates will be given per level of this factors. Change index as needed. 

slope_estimates(data = df,
                modelTMB = mod,
                num_predictor = "Weight_z",  # Name of one numeric predictor for which to estimate slopes. Change index as needed.
                fac_predictors = var_fac2[1]) # # Name of one (or two!) factor predictor(s). Slope estimates will be given per level of this factors. Change index as needed. 

slope_estimates(data = df,
                modelTMB = mod,
                num_predictor = "Station.Count_z",  # Name of one numeric predictor for which to estimate slopes. Change index as needed.
                fac_predictors = var_fac2[1]) # # Name of one (or two!) factor predictor(s). Slope estimates will be given per level of this factors. Change index as needed. 

slope_estimates(data = df,
                modelTMB = mod,
                num_predictor = "Signal.max_z",  # Name of one numeric predictor for which to estimate slopes. Change index as needed.
                fac_predictors = var_fac2[1]) # # Name of one (or two!) factor predictor(s). Slope estimates will be given per level of this factors. Change index as needed. 

#==========================================================#
# * 7.3 Estimation of differences between factor levels ---- 
#==========================================================#
# Extract pairwise comparisons among levels of specified factor predictors (means and 95% compatibility intervals).
pairwise_comparisons(data = df,
                     modelTMB = mod,
                     predictors = var_fac2[3], # Name of one (or more!) factor predictor(s).
                     component = "cond",       # Returns differences for the conditional model ("cond"), or the zero-inflation part of the model ("zi").
                     dispScale = "response",  # "response" returns absolute differences for link identity models, and ratio or odds.ratios for log or logit link models 
                                              # "link" returns estimates on the link scale
                     contrasts = "all")       # "all" returns all pairwise comparisons, "within" only the comparisons among the levels of the first factor specified within each level of the second one


#-############################################-#
# 8. Graphical display of model predictions ----
#-############################################-#
# This section produces a combined display of raw data with model-derived predictions and their 95% compatibility intervals.

#==============================================#
# * 8.1 Specify variables for the final plot ---- 
#==============================================#
# Select ONE or TWO model predictor(s) for the final plot:
plot_predictors <- c("detR", "AcSc")      # Template: Replace with the names of MAX 2 predictor variable(s).

# Select the random level of interest to display appropriate independent replicates of your raw data:
plot_random <- unique(df$Individual) #c("TT090C", "TT163C")                             # Assign NA if missing.
# plot_random <- c("rand.1", "rand.2")        # Template: Replace with the names of your random term(s).


#==================================#
# * 8.2 Derive model predictions ---- 
#==================================#
# This function generates a grid for all possible predictor combinations. It derives model predictions 
# averaged for the predictors to plot, and projects z-transformed predictors back to their raw scale.
mod_predictions <- post_predict(data = df,
                                modelTMB = mod,
                                plot_predictors = plot_predictors, # Name of predictors specified in 8.1.
                                # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                                component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").
                           
ggplot(mod_predictions) + 
  geom_ribbon(aes(x = detR, ymin = lower, ymax = upper, group = AcSc, fill = AcSc), 
              alpha = 0.5, color = NA) +
  geom_line(aes(x = detR, y = median, color = AcSc, group = AcSc), 
            alpha = 0.8, lwd = 1) +
  scale_color_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  scale_fill_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  
  ylab("position error") +
  
  ylim(0, NA) +
  
  theme_light()       
     
#plot(ggpredict(mod.s, terms = c("detR_z [all]", "pos"), bias_correction = T))                   

#============================#
# * 8.3 Summarise raw data ---- 
#============================#
# This function computes a summary of raw data at the TRUE replicate level. These values will be displayed in the final plot:
data_summary <- display_raw(data = df,
                           modelTMB = mod,
                           plot_predictors = plot_predictors, # Name of predictors specified in 8.1.
                           # plot_random = plot_random,       # (optional) - Name of the random level specified in 8.1.
                           response = var_resp,               # Name of response variable.
                           # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                           component = "all")                 # Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").



#============================#
# * 8.4 Produce final plot ---- 
#============================#
# All optional terms can be activated as desired. Default values as specified. 
final_plotting(data_summary = data_summary,       # Summary of the raw data calculated in 8.3.
               predictions = mod_predictions,    # Model predictions calculated in 8.2.
               predictors = plot_predictors,     # Name of predictors specified in 8.1.
               response = var_resp,              # Name of response variable.
               interaction.lines = T,            # (optional) - Logical (T or F). Adds interaction lines - 2 factor plots only.
               size.points = 2.5,                # (optional) - Size of the raw data points - for all plots.
               shape.points = 21,                # (optional) - Shape of the raw data points - for single-predictor plots only.
               col.fill.points = "grey50",       # (optional) - Fill colour of raw data points - for all plots.
               col.outline.points = "black",     # (optional) - contour colour of raw data points - for all plots.
               alpha.points = 0.6,               # (optional) - Transparency of raw data points - for all plots.
               jitter.points = 0.1,              # (optional) - Jitter of raw data points - for all plots.
               width.error.bars = 0.12,          # (optional) - Width of error bars - for all plots.
               size.error.bars = 0.75,           # (optional) - Size of error bars - for all plots.
               col.error.bars = "black",         # (optional) - Colour of error bars - for one predictor plots only.
               alpha.error.bars = 0.6,           # (optional) - Transparency of shaded error bars - for 2 numeric predictors plots only.
               size.means = 3.5,                 # (optional) - Size of mean points - for all factor plots.
               grouped.shape.points = NULL,      # (optional, default to shape 21) - Replace with a vector of shapes numbers of the same length as the levels of the grouping factor predictor for grouped shapes.
               grouped.errorbar.col = NULL,      # (optional, default to random colours) - Replace with a vector of colours of the same length as the levels of the grouping factor predictor for grouped error bars colours.
               grouped.fill.col = NULL,          # (optional, default to random colours) - Replace with a vector of colours of the same length as the levels of the grouping factor predictor for grouped point fill.
               grouped.linetype = NULL,          # (optional, default to lty = 1) - Replace with a vector of linetypes of the same length as the levels of the grouping factor predictor. 
               y.limits = NULL,                  # (optional, default set by ggplot) - Replace with limits for y.axis expressed as c(min,max).
               y.breaks = NULL,                  # (optional, default set by ggplot) - Replace with a vector of breaks for y.axis.
               x.lab = NULL,                     # (optional, default set by ggplot) - Replace with title for x-axis.
               y.lab = NULL,                     # (optional, default set by ggplot) - Replace with title for y-axis.
               plot.title = NULL,                # (optional) - Replace with title for the plot.
               leg.pos = "inside.top.left",      # (optional) Available for two predictors plot only. Replace with one among NULL, inside.top.right, inside.bottom.left, inside.bottom.right, outside.top, outside.left, outside.right, outside.bottom.
               leg.title = NULL,                 # (optional) Available for two predictors plots only. Replace with your legend title.
               A.ratio = NULL)                   # (optional) Define plot aspect ratio (default: 0.7).

ggplot(mod_predictions) + 
  geom_ribbon(aes(x = dens, ymin = lower, ymax = upper, group = meth, fill = meth), 
              alpha = 0.5, color = NA) +
  geom_line(aes(x = dens, y = median, color = meth, group = meth), 
            alpha = 0.8, lwd = 1) +
  scale_color_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  scale_fill_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  
  ylab("position error") +

  ylim(0, NA) +
  
  theme_light() 

ggplot(mod_predictions) + 
  geom_ribbon(aes(x = detR, ymin = lower, ymax = upper, group = pos, fill = pos), 
              alpha = 0.5, color = NA) +
  geom_line(aes(x = detR, y = median, color = pos, group = pos), 
            alpha = 0.8, lwd = 1) +
  scale_color_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  scale_fill_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  
  ylab("position error") +
  
  ylim(0, NA) +
  
  theme_light() 

ggplot(mod_predictions) + 
  geom_errorbar(aes(x = pos, ymin = lower, ymax = upper, group = meth, color = meth), 
              alpha = 0.5, color = NA) +
  geom_point(aes(x = pos, y = median, color = meth, group = meth), 
            alpha = 0.8, lwd = 1) +
  scale_color_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  scale_fill_viridis_d("method", end = 0.9, begin = 0.2, option = "rocket")+
  
  ylab("position error") +
  facet_wrap(~model) +
  
  ylim(0, NA) +
  
  theme_light() 


pred.detR <- post_predict(data = df,
                          modelTMB = mod,
                          plot_predictors = c("detR_z", "pos"), # Name of predictors specified in 8.1.
                          # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                          component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").

pred.dens <- post_predict(data = df,
                          modelTMB = mod,
                          plot_predictors = c("dens_z", "pos"), # Name of predictors specified in 8.1.
                          # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                          component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").

pred.NAnt <- post_predict(data = df,
                          modelTMB = mod,
                          plot_predictors = c("Antenna.Count_z", "pos"), # Name of predictors specified in 8.1.
                          # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                          component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").

pred.NStat <- post_predict(data = df,
                           modelTMB = mod,
                           plot_predictors = c("Station.Count_z", "pos"), # Name of predictors specified in 8.1.
                           # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                           component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").

pred.weight <- post_predict(data = df,
                            modelTMB = mod,
                            plot_predictors = c("Weight_z", "pos"), # Name of predictors specified in 8.1.
                            # offset = NA,                     # (optional) - Name of response offset variable, if present in the model. Default to NA.
                            component = "all")                 # (optional) - Computes predictions for the conditional model ("cond"), zero-inflation part of the model ("zi"), or the default both ("all").

saveRDS(mod, "./R_model/lmm_poly.RDS")

## merge predictions
colnames(pred.detR)[1] <- "pred"
pred.detR$pred.name <- "detR"
colnames(pred.dens)[1] <- "pred"
pred.dens$pred.name <- "dens"
colnames(pred.NAnt)[1] <- "pred"
pred.NAnt$pred.name <- "NAnt"
colnames(pred.NStat)[1] <- "pred"
pred.NStat$pred.name <- "NStat"
colnames(pred.weight)[1] <- "pred"
pred.weight$pred.name <- "weight"

df.pred <- rbind(pred.detR, pred.dens)
df.pred <- rbind(df.pred, pred.NAnt)
df.pred <- rbind(df.pred, pred.NStat)
df.pred <- rbind(df.pred, pred.weight)

write.csv(df.pred, "./R_model/predictions_lmm_poly.csv", row.names = F)


for(p in unique(df.pred$pred.name)) {
  
  g <- ggplot(df.pred[df.pred$pred.name == p,]) + 
    geom_ribbon(aes(x = pred, ymin = 10^lower, ymax = 10^upper, group = pos, fill = pos), 
                    alpha = 0.5, color = NA) +
    geom_line(aes(x = pred, y = 10^median, color = pos, group = pos), 
                 alpha = 0.8, lwd = 1) +
    scale_color_viridis_d("position", end = 0.9, begin = 0.2, option = "rocket")+
    scale_fill_viridis_d("position", end = 0.9, begin = 0.2, option = "rocket")+
    
    ylab("position error") +
    xlab(p) +
    
    ylim(0, NA) +
    
    theme_light()  
  
  print(g)
}

