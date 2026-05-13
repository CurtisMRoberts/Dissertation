#############################
#Chapter 5 Real-data application
#implementation on PRIME colorectal cancer trial 
#############################
rm(list=ls())


library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(geepack)
library(mice)
library(rjags)
library(coda)

#############################
#load lon PRIME dataset with VISIT labels retained
#############################

load("lon_with_visit.RData")

prime_lon <- lon

rm(lon)

#############################
#Construct scheduled-visit dataset
#############################

#The PRIME longitudinal dataset contains observed SLD records
#VISIT is retained so that scheduled visits can be indexed explicitly as j = 0, ..., T
#This allowed the missingness indicators R_ij to be constructed 

prime_lon <- prime_lon %>%
  mutate(
    VISIT = as.character(VISIT),
    y = log(SLD + 1)
  )

#Define scheduled visit grid
scheduled_visits <- c(
  "Screening",
  "Week 8",
  "Week 16",
  "Week 24",
  "Week 32",
  "Week 40",
  "Week 48",
  "Week 56",
  "Week 64",
  "Week 72",
  "Week 80",
  "Week 88",
  "Week 96"
)

#Keep responses whose VISIT label corresponds to one of the scheduled visits
#This does not exclude patients, it defines which visit labels form the scheduled grid
scheduled_lon <- prime_lon %>%
  filter(VISIT %in% scheduled_visits)

#Construct the scheduled visit grid.
#The median observed TIME for each VISIT is used as the numerical time value associated with that scheduled visit label
#this is because we have two time variables (VISIT and TIME), but TIME values arent the same for everyone so we need a single numerical value for t_j 
#to use in marginal model
visit_grid <- scheduled_lon %>%
  group_by(VISIT) %>%
  summarise(
    time = median(TIME, na.rm = TRUE),
    n_observed_values = n(),
    .groups = "drop"
  ) %>%
  arrange(time) %>%
  mutate(
    visit = row_number() - 1
  )

visit_grid #check structure

#############################
#Check duplicate observed values within subject and scheduled visit
#############################

duplicate_subject_visits <- scheduled_lon %>%
  inner_join(
    visit_grid %>% select(VISIT, visit, grid_time = time),
    by = "VISIT"
  ) %>%
  count(id, VISIT, visit) %>%
  filter(n > 1)


#n=25 duplicated observations, with there being max 2 values 
duplicate_subject_visits %>%
  summarise(
    n_duplicate_subject_visits = n(),
    max_observations_per_subject_visit = ifelse(n() > 0, max(n), 0)
  )


##############################################################
#Create one observed response per subject and scheduled visit
#################################################################

#The methods require a single response Y_ij 
#If more than one observed SLD record is present for the same subject and VISIT label, retain the record closest to the scheduled visit time
#this resolves duplicate records but does not exclude subjects

observed_scheduled <- scheduled_lon %>%
  inner_join(
    visit_grid %>% select(VISIT, visit, grid_time = time),
    by = "VISIT"
  ) %>%
  mutate(
    time_distance = abs(TIME - grid_time)
  ) %>%
  group_by(id, VISIT, visit) %>%
  arrange(time_distance, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    id,
    VISIT,
    visit,
    time = grid_time,
    TIME_observed = TIME,
    SLD,
    y,
    TRT,
    PRSURG,
    LIVERMET,
    AGE,
    SEX,
    B_WEIGHT,
    B_HEIGHT,
    B_ECOG,
    HISSUBTY,
    B_METANM,
    DIAGTYPE,
    BMMTR1
  )

#################################################################################
#Expand all subjects over the scheduled visits and define missingness indicator 
################################################################################

#Include all subjects present in the  dataset
#No restriction imposed based on observed baseline response or number of observed visits
analysis_ids <- prime_lon %>%
  distinct(id) %>%
  arrange(id)

#Subject-level covariates are taken from the dataset.
#These covariates are constant within subject 
patient_covars <- prime_lon %>%
  arrange(id, TIME) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    id,
    TRT,
    PRSURG,
    LIVERMET,
    AGE,
    SEX,
    B_WEIGHT,
    B_HEIGHT,
    B_ECOG,
    HISSUBTY,
    B_METANM,
    DIAGTYPE,
    BMMTR1
  )

#Construct the subject-by-scheduled-visit dataset
prime_full <- tidyr::expand_grid(
  id = analysis_ids$id,
  visit = visit_grid$visit
) %>%
  left_join(
    visit_grid %>% select(visit, VISIT, time),
    by = "visit"
  ) %>%
  left_join(
    observed_scheduled %>% select(id, visit, SLD, y, TIME_observed),
    by = c("id", "visit")
  ) %>%
  left_join(
    patient_covars,
    by = "id"
  ) %>%
  arrange(id, visit) %>%
  mutate(
    R = ifelse(!is.na(y), 1L, 0L)
  )

##################################################
#Construct baseline and observed-history variables
##################################################

#baseline_y is the actual Screening response. It is not replaced by the first available post-baseline response

prime_full <- prime_full %>%
  group_by(id) %>%
  arrange(visit, .by_group = TRUE) %>%
  mutate(
    baseline_y = y[visit == 0][1],
    last_obs_y = lag(zoo::na.locf(y, na.rm = FALSE)),
    last_obs_y = ifelse(visit == 0, baseline_y, last_obs_y)
  ) %>%
  ungroup()

#############################
#Checks for the constructed scheduled-visit dataset
#############################

#dimensions (n=442 indivs with 13 visits)
prime_full %>%
  summarise(
    n_subjects = n_distinct(id),
    n_visits = n_distinct(visit),
    n_rows = n()
  )

#Every subject should have one row per scheduled visit (13 rows for 13 visits as required)
prime_full %>%
  count(id) %>%
  summarise(
    min_rows_per_subject = min(n),
    median_rows_per_subject = median(n),
    max_rows_per_subject = max(n)
  )

#Overall missingness (1922 observed, 3824 missing, proportion of missingness is 0.666)
prime_full %>%
  summarise(
    n_subjects = n_distinct(id),
    n_rows = n(),
    observed = sum(R == 1),
    missing = sum(R == 0),
    prop_missing = mean(R == 0),
    min_time = min(time),
    max_time = max(time)
  )

#Missingness by scheduled visit.
prime_full %>%
  group_by(visit, VISIT, time) %>%
  summarise(
    observed = sum(R == 1),
    missing = sum(R == 0),
    prop_missing = mean(R == 0),
    .groups = "drop"
  )

#Check whether baseline response is observed, this is just a check and isnt forced to be true
prime_full %>%
  filter(visit == 0) %>%
  summarise(
    n_subjects = n_distinct(id),
    observed_baseline_y = sum(!is.na(y)),
    missing_baseline_y = sum(is.na(y)),
    prop_missing_baseline_y = mean(is.na(y))
  )

#Distribution of observed scheduled visits per subject.
obs_visit_distribution <- prime_full %>%
  group_by(id) %>%
  summarise(
    n_obs = sum(R == 1),
    has_observed_baseline = any(visit == 0 & R == 1),
    .groups = "drop"
  )

obs_visit_distribution %>%
  count(n_obs)

obs_visit_distribution %>%
  summarise(
    n_subjects = n(),
    min_observed_visits = min(n_obs),
    median_observed_visits = median(n_obs),
    max_observed_visits = max(n_obs),
    n_with_zero_observed_visits = sum(n_obs == 0),
    n_with_one_observed_visit = sum(n_obs == 1),
    prop_with_one_observed_visit = mean(n_obs == 1),
    n_missing_baseline = sum(!has_observed_baseline),
    prop_missing_baseline = mean(!has_observed_baseline)
  )



#prepare variables for modelling

prime_full <- prime_full %>%
  mutate(
    id = as.integer(id),
    R = as.integer(R),
    visit = as.integer(visit),
    time = as.numeric(time),
    
    #Make sure variables are numeric
    y = as.numeric(y),
    SLD = as.numeric(SLD),
    last_obs_y = as.numeric(last_obs_y),
    baseline_y = as.numeric(baseline_y),
    
    #Ensure baseline covariates are factors
    TRT = droplevels(factor(TRT)),
    BMMTR1 = droplevels(factor(BMMTR1)),
    LIVERMET = droplevels(factor(LIVERMET)),
    B_ECOG = droplevels(factor(B_ECOG)),
    
    #Recreate treatment indicator
    #z = 0 for FOLFOX alone
    #z = 1 for Panitumumab + FOLFOX
    z = as.numeric(TRT == "Panitumumab + FOLFOX"),
    
    #treatment and time interaction.
    time_z = time * z
  )

#Set the reference levels 
prime_full$TRT <- relevel(prime_full$TRT, ref = "FOLFOX alone")
prime_full$BMMTR1 <- relevel(prime_full$BMMTR1, ref = "Wild-type")
prime_full$LIVERMET <- relevel(prime_full$LIVERMET, ref = "N")
prime_full$B_ECOG <- relevel(prime_full$B_ECOG, ref = "Fully active")

#Check treatment balance (223 F alone, 219 P + F)
prime_full %>% 
  distinct(id, TRT) %>%
  count(TRT)

#Check KRAS status (260 wild, 182 mutant)
prime_full %>% 
  distinct(id, BMMTR1) %>%
  count(BMMTR1)

#Check liver metastasis status (N=33, Y=409 VERY UNBALANCED)
prime_full %>%
  distinct(id, LIVERMET) %>%
  count(LIVERMET)

#Check ECOG status (233 Fully active, 21 in bed less than 50% , 188 symptoms but ambulatory)
prime_full %>%
  distinct(id, B_ECOG) %>%
  count(B_ECOG)

##########################
# Structural checks
###########################

#Check dimensions of the dataset (n = 442 on 13 visits)
prime_full %>%
  summarise(
    n_subjects = n_distinct(id),
    n_visits = n_distinct(visit),
    n_rows = n()
  )

#Check that every subject has exactly one row per retained scheduled visit
prime_full %>%
  count(id) %>%
  summarise(
    min_rows_per_subject = min(n),
    median_rows_per_subject = median(n),
    max_rows_per_subject = max(n)
  )


#Check whether covariates are complete at baseline (they are)
prime_full %>%
  filter(visit == 0) %>%
  summarise(
    missing_TRT = sum(is.na(TRT)),
    missing_BMMTR1 = sum(is.na(BMMTR1)),
    missing_LIVERMET = sum(is.na(LIVERMET)),
    missing_B_ECOG = sum(is.na(B_ECOG)),
    missing_baseline_y = sum(is.na(baseline_y)),
    missing_last_obs_y = sum(is.na(last_obs_y))
  )


##########################
#Check missingness structure
###########################

#Overall missingness in the dataset
prime_full %>%
  summarise(
    n_subjects = n_distinct(id),
    n_rows = n(),
    observed = sum(R == 1),
    missing = sum(R == 0),
    prop_missing = mean(R == 0),
    min_time = min(time),
    max_time = max(time)
  )

#Missingness by scheduled visit
prime_full %>%
  group_by(visit, VISIT, time) %>%
  summarise(
    observed = sum(R == 1),
    missing = sum(R == 0),
    prop_missing = mean(R == 0),
    .groups = "drop"
  )

#Observed outcome range.
summary(prime_full$y[prime_full$R == 1])

# Quick check for the observed history variable used in WGEE and DR-GEE.
summary(prime_full$last_obs_y)
sum(is.na(prime_full$last_obs_y))


####################################################
#Save prime_full so that a privacy-preserved subset can be made
###################################################

save(prime_full, file = "prime_full.RData")


##########################
#Marginal mean model
###########################

#the common marginal mean structure used across methods
#quick note: LIVERMET is highly imbalanced in the dataset, so its coefficient should be interpreted cautiously

marginal_mean <- y ~ time * z + BMMTR1 + LIVERMET + B_ECOG


#set seed for reproducibility
set.seed(3737)


#############################
#Available-data GEE
#############################

fit_available_gee <- function(dat) {
  
  #This method only uses observed outcomes
  #does not explicitly correct for missingness
  
  fit <- geepack::geeglm(
    formula = marginal_mean,
    id = id,
    data = dat %>% filter(R == 1),
    corstr = "exchangeable"
  )
  
  beta_hat <- coef(fit)
  se_hat <- coef(summary(fit))[, "Std.err"]
  
  list(
    beta = beta_hat,
    se = se_hat,
    fit = fit,
    summary = summary(fit)
  )
}

#############################
#LOCF
#############################

fit_locf_real <- function(dat) {
  
  #LOCF replaces each missing response with the most recent previously observed response from the same individual
  #baseline is observed for all patients so this shouldnt remove any rows
  
  dat_locf <- dat %>%
    group_by(id) %>%
    arrange(visit, .by_group = TRUE) %>%
    mutate(
      y_locf = zoo::na.locf(y, na.rm = FALSE)
    ) %>%
    ungroup()
  
  #Baseline is observed, should remove no rows
  dat_locf <- dat_locf %>%
    filter(!is.na(y_locf))
  
  #Fit the same marginal mean structure, replacing y by y_locf
  fit <- geepack::geeglm(
    y_locf ~ time * z + BMMTR1 + LIVERMET + B_ECOG,
    id = id,
    data = dat_locf,
    corstr = "exchangeable"
  )
  
  beta_hat <- coef(fit)
  se_hat <- coef(summary(fit))[, "Std.err"]
  
  list(
    beta = beta_hat,
    se = se_hat,
    fit = fit,
    summary = summary(fit),
    locf_data = dat_locf
  )
}

###################
#Weighted GEE
##################

fit_wgee_real <- function(dat, min_pi = 0.05) {
  
  dat_w <- dat %>%
    arrange(id, visit)
  
  #Fit the missingness model
  #R = 1 means the response is observed at that scheduled visit
  #Baseline is excluded because it is observed for every selected patient
  #last_obs_y is included because missingness may depend on most recently previously observed response 
  
  obs_model <- glm(
    R ~ last_obs_y + z + time + BMMTR1 + LIVERMET + B_ECOG,
    data = dat_w %>% filter(visit > 0),
    family = binomial()
  )
  
  #Estimate probabilities of observation
  #and truncate very small fitted probabilities to avoid extreme weights
  
  dat_w <- dat_w %>%
    mutate(
      pi_hat_raw = ifelse(
        visit == 0,
        1,
        predict(obs_model, newdata = ., type = "response")
      ),
      pi_hat = ifelse(
        visit == 0,
        1,
        pmax(pi_hat_raw, min_pi)
      ),
      wgee_weight = ifelse(R == 1, 1 / pi_hat, NA_real_)
    )
  #Diagnostic checks for the observation probabilities and weights
  weight_diagnostics <- dat_w %>%
    filter(visit > 0) %>%
    summarise(
      min_pi_raw = min(pi_hat_raw, na.rm = TRUE),
      q1_pi_raw = quantile(pi_hat_raw, 0.25, na.rm = TRUE),
      median_pi_raw = median(pi_hat_raw, na.rm = TRUE),
      q3_pi_raw = quantile(pi_hat_raw, 0.75, na.rm = TRUE),
      max_pi_raw = max(pi_hat_raw, na.rm = TRUE),
      n_pi_below_min = sum(pi_hat_raw < min_pi, na.rm = TRUE)
    )
  
  observed_weight_diagnostics <- dat_w %>%
    filter(R == 1) %>%
    summarise(
      min_weight = min(wgee_weight, na.rm = TRUE),
      median_weight = median(wgee_weight, na.rm = TRUE),
      max_weight = max(wgee_weight, na.rm = TRUE),
      mean_weight = mean(wgee_weight, na.rm = TRUE)
    )
  
  #Fit WGEE using observed outcomes only, exchangeable correlation structure is assumed
  fit <- geepack::geeglm(
    formula = marginal_mean,
    id = id,
    data = dat_w %>% filter(R == 1),
    weights = wgee_weight,
    corstr = "exchangeable" 
  )
  
  beta_hat <- coef(fit)
  se_hat <- coef(summary(fit))[, "Std.err"]
  
  list(
    beta = beta_hat,
    se = se_hat,
    fit = fit,
    summary = summary(fit),
    obs_model = obs_model,
    wgee_data = dat_w,
    weight_diagnostics = weight_diagnostics,
    observed_weight_diagnostics = observed_weight_diagnostics
  )
}

#############################
#Multiple Imputation by Chained Equations
#############################

fit_mice_real <- function(dat, m = 20, maxit = 10, seed = 1000) {
  
  #Convert from long to wide format
  #so that each subject has one row, and repeated outcomes become y0, y1, ..., y12
  
  dat_wide <- dat %>%
    dplyr::select(id, visit, y, z, BMMTR1, LIVERMET, B_ECOG) %>%
    mutate(y_name = paste0("y", visit)) %>%
    dplyr::select(id, z, BMMTR1, LIVERMET, B_ECOG, y_name, y) %>%
    pivot_wider(
      names_from = y_name,
      values_from = y
    ) %>%
    arrange(id)
  
  #make sure covariates are factors
  dat_wide <- dat_wide %>%
    mutate(
      BMMTR1 = factor(BMMTR1),
      LIVERMET = factor(LIVERMET),
      B_ECOG = factor(B_ECOG)
    )
  
  #need to identify repeated outcome columns
  y_vars <- grep("^y[0-9]+$", names(dat_wide), value = TRUE)
  
  #Only impute outcome columns that actually contain missing values
  impute_vars <- y_vars[colSums(is.na(dat_wide[y_vars])) > 0]
  
  #Baseline outcome is observed so should not be imputed 
  impute_vars <- setdiff(impute_vars, "y0")
  
  #MICE predictor matrix
  pred <- mice::make.predictorMatrix(dat_wide)
  pred[,] <- 0
  
  #For each incomplete outcome, use:
  #treatment, baseline covariates, and all other repeated outcomes
  #This mirrors the simulation approach but may be unstable because later visits have very few observed outcomes (discuss this further)
  
  for (v in impute_vars) {
    pred[v, c("z", "BMMTR1", "LIVERMET", "B_ECOG", setdiff(y_vars, v))] <- 1
  }
  
  #Dont use id as a predictor
  pred[, "id"] <- 0
  
  #Dont impute or use id as an outcome
  pred["id",] <- 0
  
  #Gaussian imputation is common and used here because y = log(SLD + 1) is continuous
  meth <- mice::make.method(dat_wide)
  meth[] <- ""
  meth[impute_vars] <- "norm"
  
  #Run MICE
  imp <- mice::mice(
    dat_wide,
    m = m,
    maxit = maxit,
    method = meth,
    predictorMatrix = pred,
    seed = seed,
    printFlag = FALSE
  )
  
  #store estimates and covariance matrices from each imputed dataset
  beta_list <- list()
  vcov_list <- list()
  
  #map the constructed visit index to scheduled follow-up time
  time_lookup <- dat %>%
    distinct(visit, time) %>%
    arrange(visit)
  
  for (k in 1:m) {
    
    #Extract completed dataset k
    completed_wide <- mice::complete(imp, action = k)
    
    #convert completed data back to long format
    completed_long <- completed_wide %>%
      pivot_longer(
        cols = all_of(y_vars),
        names_to = "visit",
        values_to = "y"
      ) %>%
      mutate(
        visit = as.numeric(gsub("y", "", visit))
      ) %>%
      left_join(time_lookup, by = "visit") %>%
      arrange(id, visit)
    
    #Fit the marginal mean model in imputed dataset k
    fit_k <- geepack::geeglm(
      formula = marginal_mean,
      id = id,
      data = completed_long,
      corstr = "exchangeable"
    )
    
    beta_list[[k]] <- coef(fit_k)
    vcov_list[[k]] <- fit_k$geese$vbeta
  }
  
  #Pool the estimates using Rubin's rules
  beta_mat <- do.call(rbind, beta_list)
  beta_bar <- colMeans(beta_mat)
  
  W <- Reduce("+", vcov_list) / m
  B <- stats::cov(beta_mat)
  T_var <- W + (1 + 1 / m) * B
  
  se <- sqrt(diag(T_var))
  
  list(
    beta = beta_bar,
    se = se,
    mice_object = imp,
    beta_imputations = beta_mat,
    within_var = W,
    between_var = B,
    total_var = T_var,
    predictor_matrix = pred,
    method = meth
  )
}



####################
#Data Augmentation
###################

fit_da_real <- function(dat,
                        n_iter = 5000,
                        n_burn = 2000,
                        n_thin = 2,
                        n_chains = 3,
                        seed = 2000) {
  
  #Bayesian MVN model
  #Missing y values are treated as latent variables and sampled by JAGS
  #The repeated outcome vector is modelled with an exchangeable covariance structure across visits
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  dat <- dat %>%
    arrange(id, visit) %>%
    mutate(
      id_num = as.numeric(factor(id))
    )
  
  #No. of subjects and scheduled visits
  n_id <- length(unique(dat$id_num))
  T <- length(unique(dat$visit))
  
  #Convert the outcome to a subject-by-visit matrix (JAGS needs repeated outcomes in matrix form)
  dat_wide <- dat %>%
    dplyr::select(id_num, visit, y) %>%
    mutate(y_name = paste0("y", visit)) %>%
    pivot_wider(
      names_from = y_name,
      values_from = y
    ) %>%
    arrange(id_num)
  
  y_vars <- paste0("y", sort(unique(dat$visit)))
  Y <- as.matrix(dat_wide[, y_vars])
  
  #Construct design matrix for the marginal mean model
  #delete.response(terms(marginal_mean)) keeps the right-hand side only
  X_mat <- model.matrix(
    delete.response(terms(marginal_mean)),
    data = dat
  )
  
  p <- ncol(X_mat)
  beta_names <- colnames(X_mat)
  
  #Convert the design matrix into an array (subject-visit-covariate)
  X_array <- array(NA_real_, dim = c(n_id, T, p))
  
  for (q in 1:p) {
    X_array[, , q] <- matrix(
      X_mat[, q],
      nrow = n_id,
      ncol = T,
      byrow = TRUE
    )
  }
  
  #Positive-definiteness for T x T exchangeable correlation matrix needs,
  #rho > -1/(T - 1).
  rho_lower <- -1 / (T - 1) + 0.001
  rho_upper <- 0.999
  
  jags_data <- list(
    n_id = n_id,
    T = T,
    p = p,
    Y = Y,
    X = X_array,
    rho_lower = rho_lower,
    rho_upper = rho_upper
  )
  
  model_string <- "
  model {
    
    for (i in 1:n_id) {
      
      for (j in 1:T) {
        mu[i, j] <- inprod(X[i, j, ], beta[])
      }
      
      Y[i, 1] ~ dnorm(mu[i, 1], tau[1])
      resid[i, 1] <- Y[i, 1] - mu[i, 1]
      sum_resid[i, 1] <- resid[i, 1]
      
      for (j in 2:T) {
        cond_mean[i, j] <- mu[i, j] + cond_coef[j] * sum_resid[i, j - 1]
        Y[i, j] ~ dnorm(cond_mean[i, j], tau[j])
        resid[i, j] <- Y[i, j] - mu[i, j]
        sum_resid[i, j] <- sum_resid[i, j - 1] + resid[i, j]
      }
    }
    tau[1] <- 1 / pow(sigma, 2)
    
    for (j in 2:T) {
      cond_coef[j] <- rho / (1 + (j - 2) * rho)
      cond_var[j] <- pow(sigma, 2) * (1 - ((j - 1) * pow(rho, 2)) / (1 + (j - 2) * rho))
      tau[j] <- 1 / cond_var[j]
    }
    for (q in 1:p) {
      beta[q] ~ dnorm(0, 0.01)
    }
    sigma ~ dunif(0, 10)
    rho ~ dunif(rho_lower, rho_upper)
  }
  "
  
  #chain specific initial values
  inits <- function(chain = 1) {
    list(
      beta = rnorm(p, 0, 0.2),
      sigma = runif(1, 0.5, 2),
      rho = runif(1, 0.1, 0.8),
      .RNG.name = "base::Mersenne-Twister", #random number generator (Mersenne-twister is common)
      .RNG.seed = seed + chain
    )
  }
  
  init_list <- lapply(seq_len(n_chains), inits)
  
  #Compile JAGS model
  da_model <- rjags::jags.model(
    file = textConnection(model_string),
    data = jags_data,
    inits = init_list,
    n.chains = n_chains,
    n.adapt = 1000,
    quiet = TRUE
  )
  
  #burn-in period
  update(da_model, n.iter = n_burn, progress.bar = "none")
  
  #posterior sampling
  da_samples <- rjags::coda.samples(
    model = da_model,
    variable.names = c("beta", "sigma", "rho"),
    n.iter = n_iter,
    thin = n_thin,
    progress.bar = "none"
  )
  
  sample_mat <- as.matrix(da_samples)
  beta_cols <- paste0("beta[", 1:p, "]")
  beta_samples <- sample_mat[, beta_cols, drop = FALSE]
  
  beta_hat <- colMeans(beta_samples)
  se_hat <- apply(beta_samples, 2, sd)
  
  names(beta_hat) <- beta_names
  names(se_hat) <- beta_names
  
  #Gelman-Rubin convergence
  rhat <- tryCatch(
    {
      diag <- coda::gelman.diag(
        da_samples,
        autoburnin = FALSE,
        multivariate = FALSE
      )
      diag$psrf[, "Point est."]
    },
    error = function(e) {
      NULL
    }
  )
  
  list(
    beta = beta_hat,
    se = se_hat,
    da_samples = da_samples,
    beta_samples = beta_samples,
    rhat = rhat,
    model = da_model,
    Y = Y,
    X_array = X_array,
    beta_names = beta_names
  )
}

#############################
#Doubly Robust GEE
#############################

fit_drgee_real <- function(dat, min_pi = 0.05) {
  
  dat_dr <- dat %>%
    arrange(id, visit)
  
  #Missingness model same as one used in WGEE
  pi_model <- glm(
    R ~ last_obs_y + z + time + BMMTR1 + LIVERMET + B_ECOG,
    data = dat_dr %>% filter(visit > 0),
    family = binomial()
  )
  
  #Estimate and truncate observation probabilities
  dat_dr <- dat_dr %>%
    mutate(
      pi_hat_raw = ifelse(
        visit == 0,
        1,
        predict(pi_model, newdata = ., type = "response")
      ),
      
      pi_hat = ifelse(
        visit == 0,
        1,
        pmax(pi_hat_raw, min_pi)
      )
    )
  # Diagnostic checks for observation probabilities.
  pi_diagnostics <- dat_dr %>%
    filter(visit > 0) %>%
    summarise(
      min_pi_raw = min(pi_hat_raw, na.rm = TRUE),
      q1_pi_raw = quantile(pi_hat_raw, 0.25, na.rm = TRUE),
      median_pi_raw = median(pi_hat_raw, na.rm = TRUE),
      q3_pi_raw = quantile(pi_hat_raw, 0.75, na.rm = TRUE),
      max_pi_raw = max(pi_hat_raw, na.rm = TRUE),
      n_pi_below_min = sum(pi_hat_raw < min_pi, na.rm = TRUE)
    )
  
  #Outcome model used for augmentation
  #This predicts missing outcomes from observed history and baseline covariates
  outcome_model <- lm(
    y ~ time * z + BMMTR1 + LIVERMET + B_ECOG + last_obs_y,
    data = dat_dr %>% filter(R == 1, visit > 0)
  )
  
  #Construct the doubly robust pseudo-outcome
  dat_dr <- dat_dr %>%
    mutate(
      m_hat = ifelse(
        visit == 0,
        y,
        predict(outcome_model, newdata = .)
      ),
      
      y_obs_for_formula = ifelse(R == 1, y, 0),
      
      y_dr = (R / pi_hat) * y_obs_for_formula +
        (1 - R / pi_hat) * m_hat
    )
  
  # Check pseudo-outcome construction.
  dr_diagnostics <- dat_dr %>%
    summarise(
      n_missing_m_hat = sum(is.na(m_hat)),
      n_missing_y_dr = sum(is.na(y_dr)),
      min_y_dr = min(y_dr, na.rm = TRUE),
      median_y_dr = median(y_dr, na.rm = TRUE),
      max_y_dr = max(y_dr, na.rm = TRUE)
    )
  
  #Fit final marginal mean model to the pseudo-outcome
  fit <- geepack::geeglm(
    y_dr ~ time * z + BMMTR1 + LIVERMET + B_ECOG,
    id = id,
    data = dat_dr,
    corstr = "exchangeable"
  )
  
  beta_hat <- coef(fit)
  se_hat <- coef(summary(fit))[, "Std.err"]
  
  list(
    beta = beta_hat,
    se = se_hat,
    fit = fit,
    summary = summary(fit),
    pi_model = pi_model,
    outcome_model = outcome_model,
    drgee_data = dat_dr,
    pi_diagnostics = pi_diagnostics,
    dr_diagnostics = dr_diagnostics
  )
}



#############################
#Run all methods
#############################

available_result <- fit_available_gee(prime_full)

locf_result <- fit_locf_real(prime_full)

wgee_result <- fit_wgee_real(prime_full,
  min_pi = 0.05
)

mice_result <- fit_mice_real(prime_full,
  m = 20,
  maxit = 10,
  seed = 2001
)

da_result <- fit_da_real(prime_full,
  n_iter = 5000,
  n_burn = 2000,
  n_thin = 2,
  n_chains = 3,
  seed = 3001
)

drgee_result <- fit_drgee_real(prime_full,
  min_pi = 0.05
)

#############################
#Results table
#############################

make_result_df <- function(result, method) {
  
  data.frame(
    method = method,
    term = names(result$beta),
    estimate = as.numeric(result$beta),
    se = as.numeric(result$se),
    lower = as.numeric(result$beta) - 1.96 * as.numeric(result$se), #approx 95% CI's
    upper = as.numeric(result$beta) + 1.96 * as.numeric(result$se)  
  )
}

real_results <- bind_rows(
  make_result_df(available_result, "Available GEE"),
  make_result_df(locf_result, "LOCF"),
  make_result_df(wgee_result, "WGEE"),
  make_result_df(mice_result, "MICE"),
  make_result_df(da_result, "DA"),
  make_result_df(drgee_result, "DR-GEE")
)

real_results

#Rounded summary table (just used to make LaTex table easier to create)
real_results_rounded <- real_results %>%
  mutate(
    estimate = round(estimate, 3),
    se = round(se, 3),
    lower = round(lower, 3),
    upper = round(upper, 3)
  )

real_results_rounded

#################
#Diagnostics
#############

#WGEE diagnostics
wgee_result$weight_diagnostics
wgee_result$observed_weight_diagnostics
summary(wgee_result$obs_model)

#DR-GEE diagnostics
drgee_result$pi_diagnostics
drgee_result$dr_diagnostics
summary(drgee_result$pi_model)
summary(drgee_result$outcome_model)

#MICE diagnostics
#when ran MICE says "Number of logged events: "2393", mice can still produce good results but
#logged events indicate that some of the imputation models required automatic adjustment
mice_result$mice_object
mice_result$method
mice_result$predictor_matrix

#Data augmentation convergence diagnostics
da_result$rhat
summary(da_result$da_samples)
plot(da_result$da_samples)

plot(da_result$da_samples[, "rho"])
gelman.diag(da_result$da_samples[, "rho"], autoburnin = FALSE)
effectiveSize(da_result$da_samples)

#Sensitivity check for drgee and wgee with stronger truncation threshold
drgee_result_10 <- fit_drgee_real(prime_full, min_pi = 0.10)
wgee_result_10 <- fit_wgee_real(prime_full, min_pi = 0.10)

drgee_result_10$dr_diagnostics
drgee_result_10$pi_diagnostics
wgee_result_10$observed_weight_diagnostics


#################
#FIGURES
#################


#Summarise missingness at each scheduled visit.
visit_missingness <- prime_full %>%
  group_by(visit, VISIT, time) %>%
  summarise(
    n_total = n(),
    n_observed = sum(R == 1, na.rm = TRUE),
    n_missing = sum(R == 0, na.rm = TRUE),
    prop_missing = mean(R == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(visit) %>%
  mutate(
    VISIT = factor(VISIT, levels = VISIT)
  )

visit_missingness



#Missingness proportion by scheduled visit bar chart.
missingness_plot <- ggplot(
  visit_missingness,
  aes(x = VISIT, y = prop_missing)
) +
  geom_col(fill = "darkcyan", colour = "black", width = 0.7) +
  geom_text(
    aes(label = scales::percent(prop_missing, accuracy = 1)),
    vjust = -0.35,
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Scheduled visit",
    y = "Proportion missing"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_blank()
  )

missingness_plot

ggsave(
  "prime_full_missingness_prop_by_visit.png",
  plot = missingness_plot,
  width = 8,
  height = 5,
  dpi = 800
)


#Spaghetti plot

#Keep only observed responses for the plot of trajectories
spaghetti_data <- prime_full %>%
  filter(R == 1) %>%
  arrange(id, visit) %>%
  mutate(
    VISIT = factor(VISIT, levels = visit_missingness$VISIT),
    TRT = factor(TRT, levels = c("FOLFOX alone", "Panitumumab + FOLFOX"))
  )

spaghetti_plot_treatment <- ggplot(
  spaghetti_data,
  aes(x = time, y = y, group = id, colour = TRT)
) +
  #Individual trajectories
  geom_line(alpha = 0.20, linewidth = 0.65) +
  geom_point(alpha = 0.40, size = 0.45) +
  
  #treatment-specific mean trend (smooth)
  geom_smooth(
    aes(group = TRT, colour = TRT),
    method = "loess",
    se = FALSE,
    linewidth = 1.15,
    span = 0.75
  ) +
  
  scale_colour_manual(
    values = c(
      "FOLFOX alone" = "red",
      "Panitumumab + FOLFOX" = "darkcyan"
    ),
    name = "Treatment"
  ) +
  
  scale_x_continuous(
    breaks = visit_missingness$time,
    labels = as.character(visit_missingness$VISIT),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  
  labs(
    x = "Scheduled visit",
    y = expression(log(SLD + 1))
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.background = element_rect(fill = "grey70", colour = "grey40"),
    strip.text = element_text(face = "bold"),
    plot.title = element_blank()
  )

spaghetti_plot_treatment

ggsave(
  "prime_full_spaghetti_plot_treatment_colour.png",
  plot = spaghetti_plot_treatment,
  width = 8.5,
  height = 5,
  dpi = 800
)
