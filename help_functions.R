#### help functions ####

## Bayes P (P>0)
BayesP <- function(x) sum(x > 0)/length(x)

## rolling mean
rmean <- function(x, width) {
  zoo::rollapply(x,             # column to apply roll mean
                 width = width, # window size of roll mean (e.g. 10, 20, 30)
                 FUN = mean, fill = NA, na.rm = T, partial = TRUE)
}

## rolling max
rmax <- function(x, width) {
  zoo::rollapply(x,             # column to apply roll max
                 width = width, # window size of roll max (e.g. 10, 20, 30)
                 FUN = max, fill = NA, na.rm = T, partial = TRUE)
}

## function to expand to all timestamps to apply rmean and get new position error
rmeanPE <- function(data = NULL, # dataframe with signals and est. raw positions (lon, lat)
                    c.time = "X_time", # timestamp (datetime, POSIXct)
                    c.group = c("site", "tagID", "detR", "meth"), # columns used for grouping
                    c.lon = "lon", c.lat = "lat", # estimated raw positions to average
                    c.add = NULL, # additional columns you want to apply roll mean
                    GT = TRUE, # do you have ground truth data to estimate position errors?
                    c.lonT = "lon.true", c.latT = "lat.true", # groundtruth (GPS) positions, only if GT = T
                    c.PE = "PE", # position error of est. raw positions, only if GT = T
                    w.size = 30) { # window size for roll mean
  
  if(is.null(data)) stop("You need to add a dataframe.")
  
  # get dataframe
  if(GT) {
    d <- data.frame(X_time = as.POSIXct(data[, c.time]),
                    lon = data[, c.lon],
                    lat = data[, c.lat],
                    lonT = data[, c.lonT],
                    latT = data[, c.latT])    
  }
  if(!GT) {
    d <- data.frame(X_time = as.POSIXct(data[, c.time]),
                    lon = data[, c.lon],
                    lat = data[, c.lat])
  }
  
  if(!is.null(c.add)) d[, c.add] <- data[, c.add]

  d$date <- as.Date(d$X_time)
  
  d[, c.group] <- data[, c.group]
  
  df.time <- data.frame()
  
  for(i in unique(d$date)) {
    tmp <- merge(unique(d[d$date == i, c.group]),
                 data.frame(X_time = seq(from = min(d$X_time[d$date == i]), 
                                         to = max(d$X_time[d$date == i]), by = "1 sec")))
    
    df.time <- rbind(df.time, tmp)
    
  }
  
  ## merge with dfs
  d <- left_join(df.time, d, 
                     by = c(c.group, "X_time"))
  
  ## rolling mean positions
  ## here you need the full dataset (including all timestamps with NA lat and lon)
  d <- d %>% group_by(across(all_of(c.group))) %>%
    mutate(
      lon.m = rmean(lon, width = w.size),
      lat.m = rmean(lat, width = w.size),
      Weight.m = if("Weight" %in% names(d)) rmean(Weight, width = w.size) else NA_real_,
      maxSig.m = if("maxSig" %in% names(d)) rmean(maxSig, width = w.size) else NA_real_,
      maxSig.max = if("maxSig" %in% names(d)) rmax(maxSig, width = w.size) else NA_real_,
      Ac.m = if("Ac" %in% names(d)) rmean(Ac, width = w.size) else NA_real_,
      Sc.m = if("Sc" %in% names(d)) rmean(Sc, width = w.size) else NA_real_,
    ) %>%
    ungroup()   
  
  ## Position Error (PE) based on mean positions
  ## here you need only data != na in lon lat (to compute distances)
  d <- d[!is.na(d$lon),]
  
  if(GT) {
    d <- d %>% rowwise %>%
      mutate(PE = distm(x = c(lon.m, lat.m),
                        y = c(lonT, latT))) %>%
      ungroup()    
  }
 
  ## remove or rename duplicate columns to generate identical df to data (but with mean lon lat and new PE)
  d[, c.lon]<- d$lon.m
  d[, c.lat]<- d$lat.m
  colnames(d)[colnames(d) == "X_time"] <- c.time
  colnames(d)[colnames(d) == "PE"] <- c.PE
  
  if(GT) {
    data <- dplyr::select(data, -all_of(c(c.lon, c.lat, c.PE)))
    d <- dplyr::select(d, -c("lon.m", "lat.m", "lonT", "latT", "date"))
  }
  if(!GT) {
    data <- dplyr::select(data, -all_of(c(c.lon, c.lat)))
    d <- dplyr::select(d, -c("lon.m", "lat.m", "date"))
  }
  
  data <- left_join(data, d, by = c(c.group, c.time))
  
  return(data)
}

## simulate raw data to get quantiles for newdat (e.g., median, 65%)
quant.rlnorm <- function(m, # matrix with mean
                         sd, # matrix with sd
                         probs = 0.5, # quantile to calculate (0.5, 0.65, 0.75, ...)
                         nsim = 1000, # number of simulations per dataset for rlnorm
                         SPEED = F # do you want to speed up the process using parallelization (TRUE)?
) {
  
  if(length(probs) > 1) {
    warning("There is more than one value for 'probs'. The first one will be selected for further calculations.")
    probs <- probs[1]
  }
  
  # v <- sd^2
  # phi <- sqrt(v+m^2)
  # mu <- log(m^2/phi)
  sigma <- sqrt(log(1+sd^2/m^2))
  mu <- log(m)-0.5*sigma^2 # identisch zu oberem mu
  
  quant.rlnorm2 <- function(mus) {
    set.seed(3171)
    quantile(rlnorm(nsim, mus[1], mus[2]), probs = probs)  # change to rnorm, rpois, etc. if needed
  }
  
  if(SPEED) {
    require(future.apply)
    plan(multicore, workers = parallel::detectCores())
    q <- matrix(future_apply(cbind(as.vector(mu), as.vector(sigma)), 1, 
                             quant.rlnorm2, future.seed = 3171), 
                nrow = nrow(mu))
    
    plan(sequential)
    
  }
  
  if(!SPEED) {
    q <- matrix(apply(cbind(as.vector(mu), as.vector(sigma)), 1, 
                      quant.rlnorm2), 
                nrow = nrow(mu))
    
  } 
  
  return(q)
  
}

## get model predictions, adapted from 'A versatile workflow for linear modelling in R' by
#  Matteo Santon, Fraenzi Korner-Nievergelt, Nico Michiels, Nils Anthes (2023)
pred.glmmTMB <- function(data, mod, newdat, sims = F, nsim = 10000, DISP = F) { 
  
  # # Make sure table comes as dataframe, and all character vectors are factors (sometimes this causes issues):
  # data <- as.data.frame(data)
  # data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], 
  #                                            as.factor)
  
  set.seed(3171) # reproducible results
  
  formula <- paste(gsub(".~", "", formula(mod, fixed.only = T))[3], collapse = "")
  
  temp <- unlist(str_split(formula, "\\+"))
  temp <- unlist(str_split(temp, "\\*"))         # solves interaction if expressed with "\\*"
  temp <- str_replace_all(temp, fixed(" "), "")  # remove ""
  temp <- gsub("s\\(|\\)", "", temp)  # remove "s()
  
  remove_me <- temp[str_detect(temp, pattern = (":"))] # extract interactions 
  
  if(length(remove_me) > 0) {
    fixed.pred <- temp[-match(remove_me, temp)]
  } else {
    fixed.pred <- temp
  }
  
  poly_form <- fixed.pred[str_detect(fixed.pred, pattern = "poly")] # find poly formula if present
  
  if(length(poly_form) > 0) {
    
    poly_col <- gsub(" ", "", gsub(",", "", gsub("poly\\(|\\)|x\\=|[[:digit:]])||degree = [[:digit:]],|,raw=TRUE|raw=T|raw=FALSE|raw=F,|, |degree=[[:digit:]],|, |[[:digit:]]", "", poly_form)))
    
    #OLD version: gsub(" ", "", gsub("poly\\(|\\)|x\\=|[[:digit:]])||degree = [[:digit:]],|,|degree=[[:digit:]],|, |[[:digit:]],|, raw=TRUE|raw=T|raw=FALSE|raw=F", "", poly_form))
    
    for (i in 1:length(poly_form)) {
      fixed.pred <- str_replace_all(fixed.pred, fixed(poly_form[i]), poly_col[i])
    }
  } # change name of the formula to real column in list of fixed predictors
  
  # Now we generate model predictions from simulated data for data ranges defined in the newdat
  # The procedure follows the one presented in GLMMTMB paper from Brooks 2017
  
  ## get simulations
  if(DISP) {
    bsim_disp <- mvrnorm(nsim, 
                         mu = fixef(mod)$disp, 
                         Sigma = vcov(mod)$disp) # equivalent to sim from package arm.
    formula <- formula(mod$modelInfo$allForm$dispformula)
    
    Xmat_disp <- model.matrix(formula, data = newdat)
    fitmatrix_disp <- mod$modelInfo$family$linkinv(Xmat_disp %*% t(bsim_disp)) 
    pred_disp  <- mod$modelInfo$family$linkinv(Xmat_disp %*% fixef(mod)$disp)
  }
  
  bsim_cond <- mvrnorm(nsim, 
                       mu = fixef(mod)$cond, 
                       Sigma = vcov(mod)$cond) # equivalent to sim from package arm.
  formula <- eval(str2expression(paste(gsub(".~","", formula(mod, fixed.only = T))[c(1,3)], collapse = "")))
  
  
  Xmat_cond <- model.matrix(formula, data = newdat)
  fitmatrix_all <- mod$modelInfo$family$linkinv(Xmat_cond %*% t(bsim_cond)) 
  pred_mod  <- mod$modelInfo$family$linkinv(Xmat_cond %*% fixef(mod)$cond)
  
  newdat$median <- apply(fitmatrix_all, 1, quantile, prob = 0.5)  # get median
  newdat$pred_mod <-  pred_mod
  newdat$lwr <- apply(fitmatrix_all, 1, quantile, prob = 0.025) # get lower 95% CI
  newdat$lwr25 <- apply(fitmatrix_all, 1, quantile, prob = 0.25) # get lower 50% CI
  newdat$upr75 <- apply(fitmatrix_all, 1, quantile, prob = 0.75) # get upper 50% CI
  newdat$upr <- apply(fitmatrix_all, 1, quantile, prob = 0.975) # get upper 95% CI
  newdat$SD <- apply(fitmatrix_all, 1, sd) # get SD
  
  if(DISP) {
    newdat$medianD <- apply(fitmatrix_disp, 1, quantile, prob = 0.5)  # get median
    newdat$pred_modD <-  pred_disp
    newdat$lwrD <- apply(fitmatrix_disp, 1, quantile, prob = 0.025) # get lower 95% CI
    newdat$lwr25D <- apply(fitmatrix_disp, 1, quantile, prob = 0.25) # get lower 50% CI
    newdat$upr75D <- apply(fitmatrix_disp, 1, quantile, prob = 0.75) # get upper 50% CI
    newdat$uprD <- apply(fitmatrix_disp, 1, quantile, prob = 0.975) # get upper 95% CI
    newdat$SDD <- apply(fitmatrix_disp, 1, sd) # get SD
  }
  
  # Project z-transformed predictors to raw scale:
  if(length(grep("_z", names(newdat))) > 0) {
    z_predictors <- names(newdat)[grep("_z", names(newdat))]
    predictors <- str_replace_all(z_predictors[grep("_z", z_predictors)], fixed("_z"), "")
    
    temp2 <- NULL
    for(i in 1:(length(predictors))) {
      temp1 <- back_z_transform(data[,predictors[i]], newdat[,z_predictors[i]]) 
      temp2 <- as.data.frame(cbind(temp2, temp1))
    }
    
    newdat[, z_predictors] <- temp2
    names(newdat) <- str_replace_all(names(newdat), fixed("_z"), "")
  }

  ## save simulations if sims == TRUE 
  if(sims & DISP) mod.sims <<- list(fitmatrix_all, fitmatrix_disp)
  if(sims & !DISP) mod.sims <<- list(fitmatrix_all)
  
  return(newdat)
  
}

## scale numeric predictors to mean = 0 and SD = 1, from Santon et al. 2023
z_transform <- function(data, predictors) {
  
  data <- as.data.frame(data)
  
  predictors <- na.omit(predictors)
  
  if(length(predictors) == 0) stop("ERROR: No NUMERIC predictor has been specified. Skip this step if none is available.")
  
  if(sum(sapply(data[predictors], FUN = is.numeric)) != length(predictors)) 
  { stop("ERROR: You either have no NUMERIC predictors (skip this step), or some specified predictors are not numeric (specify numeric predictors only).")} else {
    
    num <- colnames(dplyr::select_if(data[predictors], is.numeric))
    
    temp_list <- lapply(as.list(data[num]), scale)
    names(temp_list) <- paste(num,"_z", sep = "") 
    return(temp_list)
  }
}

## rescale numeric predictors, from Santon et al. 2023
back_z_transform <- function(predictor, z_predictor) {
  (z_predictor * sd(predictor)) + mean(predictor)
}

## posterior predictive checks: mean, adapted from Santon et al 2023

mean.sim <- function(data, mod, response, predictor, n.sim) {
  
  data <- as.data.frame(data)
  
  sims <- simulate(mod, nsim = n.sim, seed = 9) 
  
  if(!missing(predictor) && !is.null(predictor) && length(na.omit(predictor)) > 0) { 
    
    if(length(na.omit(predictor)) > 1) {warning(call. = F, "Multiple predictors detected: the first has been taken by default. Explicitly specify others if preferred.")}
    
    if(length(predictor) > 1) {
      predictor <- na.omit(predictor)[1]
    }
    
    if(sum(sapply(data[predictor], FUN = is.numeric)) == length(predictor))  {stop("ERROR: The predictor must be a factor.")} else {
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        list_sims <- lapply(sims, function(x) {
          cbind("estimate" = x[,1]/(x[,1] + x[,2]), data %>% dplyr::select(!!sym(predictor)))
        }
        )}
      else {
        list_sims <- lapply(sims, function(x) { # lapply applies function to each vector and returns it as a list.
          cbind("estimate" = x, data %>% dplyr::select(!!sym(predictor)))
        }
        )}
      
      sims_summary <- lapply(list_sims, function(x) {
        as.data.frame(x %>% group_by(!!sym(predictor)) %>% 
                        summarize(mean = round(mean(estimate), digits = 5)))
      }
      )
      all_sims <- do.call(rbind, sims_summary)
      temp <- data %>% group_by(!!sym(predictor)) %>% 
        summarize(mean = mean(!!sym(response)))
      
      # The following checks the quantile of observed parameters against simulated ones.
      # Depending on quantiles, observed vlines are coloured blue, orange, or red.
      
      grouped_sims <- cbind(all_sims, temp[2])
      names(grouped_sims) <- c(predictor, "mean", "obs_mean")
      
      grouped_sims <- grouped_sims %>% group_by(!!sym(predictor))
      
    }}  else { # not grouped by a predictor
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        list_sims <- lapply(sims, function(x) {
          "estimate" = x[,1]/(x[,1] + x[,2])
        }
        )}
      else {
        list_sims <- lapply(sims, function(x) {
          "estimate" = x
        }
        )}
      
      sims_summary <- lapply(list_sims, function(x) {
        as.data.frame(list(mean = round(mean(x), digits = 5)))
      })
      all_sims <- do.call(rbind, sims_summary) # do.call constructs and executes a function call from a name or 
      temp <- data %>% summarize(mean = mean(!!sym(response)))
      
      # The following checks the quantile of observed parameters against simulated ones.
      # Depending on quantiles, observed vlines are coloured blue, orange, or red.
      
      grouped_sims <- cbind(all_sims, temp)
      names(grouped_sims) <- c("mean","obs_mean")
    }
  
  prop <- as.data.frame(grouped_sims %>% summarize(extreme_mean = sum(mean > obs_mean)/n.sim, equal_mean = sum(mean == obs_mean)/n.sim))
  
  col.vline_mean <- ifelse(between(prop$extreme_mean, 0.05, 0.95) | prop$equal_mean > 0.1, "blue",  
                           ifelse(between(prop$extreme_mean, 0.005, 0.995) | prop$equal_mean > 0.01, "orange", "red"))
  
  plot1 <- ggplot() + 
    labs(x = "mean: simulated (bars) and observed (line)", y = "Frequency", title = "") +
    geom_histogram(data = all_sims, aes(x = mean), col = "white", fill = "grey70", bins = 30) + 
    geom_vline(data = temp, aes(xintercept = mean), col = col.vline_mean, linewidth = 1.5) +
    plot0 +
    if(!missing(predictor) && !is.null(predictor)  && length(na.omit(predictor)) > 0) {facet_wrap(as.formula(paste("~", predictor)), scales = "free")}
  
  if(sum(col.vline_mean == "blue")  == length(col.vline_mean)) {
    message("OPTIMAL FIT: All observed means fall within central 90% of simulated data. No issues.")} else {
      if(sum(col.vline_mean == "orange") == length(col.vline_mean)) {
        warning(call. = F, "SUB-OPTIMAL FIT: All observed means (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")} else {
          if(sum(col.vline_mean == "red") == length(col.vline_mean)) {
            warning(call. = F, "POOR FIT: All observed means (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
              if(sum(col.vline_mean == "blue") + sum(col.vline_mean == "orange") == length(col.vline_mean)) {
                warning(call. = F, "SUB-OPTIMAL FIT: Some observed means (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")} else {
                  if(sum(col.vline_mean == "blue") + sum(col.vline_mean == "red") == length(col.vline_mean)) {
                    warning(call. = F, "POOR FIT: Some observed means (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
                      if(sum(col.vline_mean == "orange") + sum(col.vline_mean == "red") == length(col.vline_mean)) {
                        warning(call. = F, "SUB-OPTIMAL FIT: Some observed means (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")
                        warning(call. = F, "POOR FIT: Some observed means (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
                          if(sum(col.vline_mean == "blue") + sum(col.vline_mean == "orange") + sum(col.vline_mean == "red") == length(col.vline_mean)) {   
                            warning(call. = F, "SUB-OPTIMAL FIT: Some observed means (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")
                            warning(call. = F, "POOR FIT: Some observed means (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")
                          }}}}}}}
  print(plot1)
  
}

## posterior predictive checks: median, adapted from Santon et al 2023

median.sim <- function(data, mod, response, predictor, n.sim) {
  
  data <- as.data.frame(data)
  
  sims <- simulate(mod, nsim = n.sim, seed = 9) 
  
  if(!missing(predictor) && !is.null(predictor) && length(na.omit(predictor)) > 0) { 
    
    if(length(na.omit(predictor)) > 1) {warning(call. = F, "Multiple predictors detected: the first has been taken by default. Explicitly specify others if preferred.")}
    
    if(length(predictor) > 1) {
      predictor <- na.omit(predictor)[1]
    }
    
    if(sum(sapply(data[predictor], FUN = is.numeric)) == length(predictor))  {stop("ERROR: The predictor must be a factor.")} else {
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        list_sims <- lapply(sims, function(x) {
          cbind("estimate" = x[,1]/(x[,1] + x[,2]), data %>% dplyr::select(!!sym(predictor)))
        }
        )}
      else {
        list_sims <- lapply(sims, function(x) { # lapply applies function to each vector and returns it as a list.
          cbind("estimate" = x, data %>% dplyr::select(!!sym(predictor)))
        }
        )}
      
      sims_summary <- lapply(list_sims, function(x) {
        as.data.frame(x %>% group_by(!!sym(predictor)) %>% 
                        summarize(median = round(median(estimate), digits = 5)))
      }
      )
      all_sims <- do.call(rbind, sims_summary)
      temp <- data %>% group_by(!!sym(predictor)) %>% 
        summarize(median = median(!!sym(response)))
      
      # The following checks the quantile of observed parameters against simulated ones.
      # Depending on quantiles, observed vlines are coloured blue, orange, or red.
      
      grouped_sims <- cbind(all_sims, temp[2])
      names(grouped_sims) <- c(predictor, "median", "obs_median")
      
      grouped_sims <- grouped_sims %>% group_by(!!sym(predictor))
      
    }}  else { # not grouped by a predictor
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        list_sims <- lapply(sims, function(x) {
          "estimate" = x[,1]/(x[,1] + x[,2])
        }
        )}
      else {
        list_sims <- lapply(sims, function(x) {
          "estimate" = x
        }
        )}
      
      sims_summary <- lapply(list_sims, function(x) {
        as.data.frame(list(median = round(median(x), digits = 5)))
      })
      all_sims <- do.call(rbind, sims_summary) # do.call constructs and executes a function call from a name or 
      temp <- data %>% summarize(median = median(!!sym(response)))
      
      # The following checks the quantile of observed parameters against simulated ones.
      # Depending on quantiles, observed vlines are coloured blue, orange, or red.
      
      grouped_sims <- cbind(all_sims, temp)
      names(grouped_sims) <- c("median","obs_median")
    }
  
  prop <- as.data.frame(grouped_sims %>% summarize(extreme_median = sum(median > obs_median)/n.sim, equal_median = sum(median == obs_median)/n.sim))
  
  col.vline_median <- ifelse(between(prop$extreme_median, 0.05, 0.95) | prop$equal_median > 0.1, "blue",  
                             ifelse(between(prop$extreme_median, 0.005, 0.995) | prop$equal_median > 0.01, "orange", "red"))
  
  plot1 <- ggplot() + 
    labs(x = "median: simulated (bars) and observed (line)", y = "Frequency", title = "") +
    geom_histogram(data = all_sims, aes(x = median), col = "white", fill = "grey70", bins = 30) + 
    geom_vline(data = temp, aes(xintercept = median), col = col.vline_median, linewidth = 1.5) +
    plot0 +
    if(!missing(predictor) && !is.null(predictor)  && length(na.omit(predictor)) > 0) {facet_wrap(as.formula(paste("~", predictor)), scales = "free")}
  
  if(sum(col.vline_median == "blue")  == length(col.vline_median)) {
    message("OPTIMAL FIT: All observed medians fall within central 90% of simulated data. No issues.")} else {
      if(sum(col.vline_median == "orange") == length(col.vline_median)) {
        warning(call. = F, "SUB-OPTIMAL FIT: All observed medians (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")} else {
          if(sum(col.vline_median == "red") == length(col.vline_median)) {
            warning(call. = F, "POOR FIT: All observed medians (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
              if(sum(col.vline_median == "blue") + sum(col.vline_median == "orange") == length(col.vline_median)) {
                warning(call. = F, "SUB-OPTIMAL FIT: Some observed medians (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")} else {
                  if(sum(col.vline_median == "blue") + sum(col.vline_median == "red") == length(col.vline_median)) {
                    warning(call. = F, "POOR FIT: Some observed medians (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
                      if(sum(col.vline_median == "orange") + sum(col.vline_median == "red") == length(col.vline_median)) {
                        warning(call. = F, "SUB-OPTIMAL FIT: Some observed medians (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")
                        warning(call. = F, "POOR FIT: Some observed medians (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
                          if(sum(col.vline_median == "blue") + sum(col.vline_median == "orange") + sum(col.vline_median == "red") == length(col.vline_median)) {   
                            warning(call. = F, "SUB-OPTIMAL FIT: Some observed medians (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")
                            warning(call. = F, "POOR FIT: Some observed medians (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")
                          }}}}}}}
  print(plot1)
  
}

## posterior predictive checks: SD, adapted from Santon et al 2023

sd.sim <- function(data, mod, response, predictor, n.sim) {
  
  data <- as.data.frame(data)
  
  sims <- simulate(mod, nsim = n.sim, seed = 9) 
  
  if(!missing(predictor) && !is.null(predictor) && length(na.omit(predictor)) > 0) { 
    
    if(length(na.omit(predictor)) > 1) {warning(call. = F, "Multiple predictors detected: the first has been taken by default. Explicitly specify others if preferred.")}
    
    if(length(predictor) > 1) {
      predictor <- na.omit(predictor)[1]
    }
    
    if(sum(sapply(data[predictor], FUN = is.numeric)) == length(predictor))  {stop("ERROR: The predictor must be a factor.")} else {
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        list_sims <- lapply(sims, function(x) {
          cbind("estimate" = x[,1]/(x[,1] + x[,2]), data %>% dplyr::select(!!sym(predictor)))
        }
        )}
      else {
        list_sims <- lapply(sims, function(x) { # lapply applies function to each vector and returns it as a list.
          cbind("estimate" = x, data %>% dplyr::select(!!sym(predictor)))
        }
        )}
      
      sims_summary <- lapply(list_sims, function(x) {
        as.data.frame(x %>% group_by(!!sym(predictor)) %>% 
                        summarize(sd = round(sd(estimate), digits = 5)))
      }
      )
      all_sims <- do.call(rbind, sims_summary)
      temp <- data %>% group_by(!!sym(predictor)) %>% 
        summarize(sd = sd(!!sym(response)))
      
      # The following checks the quantile of observed parameters against simulated ones.
      # Depending on quantiles, observed vlines are coloured blue, orange, or red.
      
      grouped_sims <- cbind(all_sims, temp[2])
      names(grouped_sims) <- c(predictor, "sd", "obs_sd")
      
      grouped_sims <- grouped_sims %>% group_by(!!sym(predictor))
      
    }}  else { # not grouped by a predictor
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        list_sims <- lapply(sims, function(x) {
          "estimate" = x[,1]/(x[,1] + x[,2])
        }
        )}
      else {
        list_sims <- lapply(sims, function(x) {
          "estimate" = x
        }
        )}
      
      sims_summary <- lapply(list_sims, function(x) {
        as.data.frame(list(sd = round(sd(x), digits = 5)))
      })
      all_sims <- do.call(rbind, sims_summary) # do.call constructs and executes a function call from a name or 
      temp <- data %>% summarize(sd = sd(!!sym(response)))
      
      # The following checks the quantile of observed parameters against simulated ones.
      # Depending on quantiles, observed vlines are coloured blue, orange, or red.
      
      grouped_sims <- cbind(all_sims, temp)
      names(grouped_sims) <- c("sd","obs_sd")
    }
  
  prop <- as.data.frame(grouped_sims %>% summarize(extreme_sd = sum(sd > obs_sd)/n.sim, equal_sd = sum(sd == obs_sd)/n.sim))
  
  col.vline_sd <- ifelse(between(prop$extreme_sd, 0.05, 0.95) | prop$equal_sd > 0.1, "blue",  
                         ifelse(between(prop$extreme_sd, 0.005, 0.995) | prop$equal_sd > 0.01, "orange", "red"))
  
  plot1 <- ggplot() + 
    labs(x = "sd: simulated (bars) and observed (line)", y = "Frequency", title = "") +
    geom_histogram(data = all_sims, aes(x = sd), col = "white", fill = "grey70", bins = 30) + 
    geom_vline(data = temp, aes(xintercept = sd), col = col.vline_sd, linewidth = 1.5) +
    plot0 +
    if(!missing(predictor) && !is.null(predictor)  && length(na.omit(predictor)) > 0) {facet_wrap(as.formula(paste("~", predictor)), scales = "free")}
  
  if(sum(col.vline_sd == "blue")  == length(col.vline_sd)) {
    message("OPTIMAL FIT: All observed sds fall within central 90% of simulated data. No issues.")} else {
      if(sum(col.vline_sd == "orange") == length(col.vline_sd)) {
        warning(call. = F, "SUB-OPTIMAL FIT: All observed sds (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")} else {
          if(sum(col.vline_sd == "red") == length(col.vline_sd)) {
            warning(call. = F, "POOR FIT: All observed sds (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
              if(sum(col.vline_sd == "blue") + sum(col.vline_sd == "orange") == length(col.vline_sd)) {
                warning(call. = F, "SUB-OPTIMAL FIT: Some observed sds (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")} else {
                  if(sum(col.vline_sd == "blue") + sum(col.vline_sd == "red") == length(col.vline_sd)) {
                    warning(call. = F, "POOR FIT: Some observed sds (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
                      if(sum(col.vline_sd == "orange") + sum(col.vline_sd == "red") == length(col.vline_sd)) {
                        warning(call. = F, "SUB-OPTIMAL FIT: Some observed sds (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")
                        warning(call. = F, "POOR FIT: Some observed sds (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")} else {
                          if(sum(col.vline_sd == "blue") + sum(col.vline_sd == "orange") + sum(col.vline_sd == "red") == length(col.vline_sd)) {   
                            warning(call. = F, "SUB-OPTIMAL FIT: Some observed sds (orange) fall marginal within central 99% of simulated data. Check if better model fit can be achieved.")
                            warning(call. = F, "POOR FIT: Some observed sds (red) are outside the central 99% of simulated data. Your model likely provides insufficient fit to your data.")
                          }}}}}}}
  print(plot1)
  
}

## posterior predictive checks: data distribution, adapted from Santon et al 2023
ppcheck_fun <- function(data, mod, response, predictor, n.sim) {
  
  data <- as.data.frame(data)
  
  sims <- simulate(mod, nsim = n.sim, seed = 9) # simulate data from the model
  
  if(!missing(predictor) && !is.null(predictor) && length(na.omit(predictor)) > 0) { 
    
    if(length(na.omit(predictor)) > 1) {warning(call. = F, "Multiple predictors detected: the first has been taken by default. Explicitly specify others if preferred.")}
    
    if(length(predictor) > 1) {
      predictor <- na.omit(predictor)[1]
    }
    
    if(sum(sapply(data[predictor], FUN = is.numeric)) == length(predictor)) {stop("ERROR: The predictor must be a factor")} else {
      
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        simulations <- lapply(sims, function(x) {
          "value" = x[,1]/(x[,1] + x[,2])})}
      else {
        simulations <- lapply(sims, function(x) {
          "value" = x})
      }
      simulations <- cbind(reshape::melt(simulations), data %>% dplyr::select(!!sym(predictor)), row.names = NULL)
    }} else {
      if (mod$modelInfo$family[1] == "binomial" | mod$modelInfo$family[1] == "betabinomial") {
        simulations <- lapply(sims, function(x) {
          "value" = x[,1]/(x[,1] + x[,2])})
      } else {
        simulations <- lapply(sims, function(x) {
          "value" = x})
      }
      simulations <- reshape::melt(simulations)
    }
  
  plot.output <- ggplot() + 
    geom_density(data = simulations, 
                 aes(x = value, group = L1), 
                 col = "grey60") + # is geom_density what we are looking for?
    geom_density(data = data, 
                 aes(x = !!sym(response)), 
                 col = "blue", linewidth = 1.5) +
    # scale_x_continuous(limits = c(0, quantile(simulations$value, prob = 0.999))) +
    labs(x = paste(response, ": simulated (grey) and observed (blue)", sep = ""), 
         y = "Density", 
         title = "") + 
    plot0
  
  if(!missing(predictor) && !is.null(predictor) && length(na.omit(predictor)) > 0) {
    plot.output <- plot.output + facet_wrap(as.formula(paste("~", predictor)), scales = "free")
  }
  print(plot.output)
}

# Theme with a legend from Santon et al. 2023
plot0.1 <- theme(axis.title.x = 
                   element_text(
                     margin = margin(t = 10),
                     color = "Black", 
                     size = 17), 
                 axis.title.y = 
                   element_text(
                     margin = margin(r = 20),
                     color = "Black", 
                     size = 17), 
                 axis.text = 
                   element_text(
                     color = "Black",
                     size = 13), 
                 axis.line = 
                   element_line(
                     color = "Black",
                     linewidth = 0.5), 
                 axis.ticks = element_line(color = "Black"),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.minor = element_line(color = "White"),
                 plot.title = 
                   element_text(
                     size = 17,
                     face = "bold",
                     hjust = 0),
                 plot.subtitle = 
                   element_text(
                     size = 15,
                     hjust = 0.5),
                 strip.text.x = 
                   element_text(
                     size = 15,
                     color = "black", 
                     face = "bold"),
                 strip.background = 
                   element_rect(
                     color = "white",
                     fill = "white", 
                     linewidth = 1, 
                     linetype = "solid"),
                 legend.background = element_blank(), 
                 legend.title = 
                   element_text(
                     colour = "black",
                     size = 13, 
                     face = "bold"), 
                 legend.text = 
                   element_text(
                     colour = "black",
                     size = 12),
                 legend.key.width = unit(1.5, "lines"), 
                 legend.key = element_rect(fill = NA))

# Same theme without a legend:
plot0 <- plot0.1 + theme(legend.position = "none")

## get coefficient estimates, from Santon et al. 2023

comp_int <- function(mod, ci_range, effects, component) {
  
  if(missing(ci_range)) ci_range = 0.95
  if(missing(effects)) effects = "all"
  if(missing(component)) component = "all"
  
  if(effects == "random") {
    message("NOTE: you selected to display random effects, only. This 
will generate an error message in case your model does not contain any 
random effect.")}
  # We use the parameters function (from package 'parameters') to extract model coefficients:
  temp <- parameters(mod,
                     ci = ci_range,
                     iterations = 10000,
                     effects = effects,
                     component = component)
  
  # We simplify the output to effect estimates and their SE / CI, but remove null hypothesis test info:
  if(effects == "all") {
    output.cols <- c("Parameter", "Coefficient", "SE", "CI_low", "CI_high", "Effects", "Group", "Component")     # c(1:3,5:6,10:12)
    if(length(temp) < 12) {  # length can be 10 or 11 in the absence of random component.
      output.cols <- c("Parameter", "Coefficient", "SE", "CI_low", "CI_high", "Effects")  # c(1:3,5:6,10)
      message("HINT: Model contains no random component")
    }
  }
  
  if(effects == "random") {
    output.cols <- c("Parameter", "Coefficient", "SE", "CI_low", "CI_high", "Effects", "Group", "Component")   # c(1:3,5:9)
  }
  if(effects == "fixed") {
    output.cols <- c("Parameter", "Coefficient", "SE", "CI_low", "CI_high", "Effects")  # c(1:3,5:6, 10)
  }
  print(temp[, output.cols])
  
}


## rolling mean
rmean <- function(x, width) {
  zoo::rollapply(x, width = width, 
                 FUN = mean, fill = NA, na.rm = T, partial = TRUE)
}
