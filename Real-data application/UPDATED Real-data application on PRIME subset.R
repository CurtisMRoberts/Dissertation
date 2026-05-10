#############################
#Chapter 5 Real-data application
#implementation on PRIME 50 patient subset
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
#load analysis dataset
#############################

prime_50 <- read.csv("prime_50_scheduled_analysis.csv")


#prepare variables for modelling

prime_50 <- prime_50 %>%
  mutate(
    #make sure the response indicator is numeric/integer
    R = as.integer(R),
    
    #make sure visit and time are numeric
    visit = as.integer(visit),
    time = as.numeric(time),
    
    #ensure categorical baseline covariates are factors
    TRT = droplevels(factor(TRT)),
    BMMTR1 = droplevels(factor(BMMTR1)),
    LIVERMET = droplevels(factor(LIVERMET)),
    B_ECOG = droplevels(factor(B_ECOG)),
    
    #recreate treatment indicator
    #z=0 for FOLFOX alone
    #z=1 for Panitumumab + FOLFOX
    z = as.numeric(TRT == "Panitumumab + FOLFOX"),
    
    #treatment and time interaction
    time_z = time * z
  )

#Set reference levels 
prime_50$TRT <- relevel(prime_50$TRT, ref = "FOLFOX alone")
prime_50$BMMTR1 <- relevel(prime_50$BMMTR1, ref = "Wild-type")


#Check treatment balance (equal)
prime_50 %>% 
  distinct(id, TRT) %>%
  count(TRT)

#Check KRAS status balance (wild=26, mutant=24)
prime_50 %>% 
  distinct(id, BMMTR1) %>%
  count(BMMTR1)

#Check liver metastasis status balance (N=3, Y=47 (HIGHLY UNBALANCED))
prime_50 %>%
  distinct(id, LIVERMET) %>%
  count(LIVERMET)

#check ECOG status (Fully active = 22, symptoms but ambulatory = 28)
prime_50 %>%
  distinct(id, B_ECOG) %>%
  count(B_ECOG)


##########################
#Check missingness structure
###########################

#Overall missingness (242 obs, 408 miss, 0.628 prop of miss)
prime_50 %>%
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
prime_50 %>%
  group_by(visit, VISIT, time) %>%
  summarise(
    observed = sum(R == 1),
    missing = sum(R == 0),
    prop_missing = mean(R == 0),
    .groups = "drop"
  )

#observed outcome range.
summary(prime_50$y[prime_50$R == 1])

#quick check for the observed history variable used in WGEE and DR-GEE 
summary(prime_50$last_obs_y)
sum(is.na(prime_50$last_obs_y))


##########################
#Marginal mean model
###########################

#the common marginal mean structure used across methods
#quick note: LIVERMET is highly imbalanced in this subset, so its coefficient should be interpreted cautiously

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
  
  dat_locf <- dat %>%
    group_by(id) %>%
    arrange(visit, .by_group = TRUE) %>%
    mutate(
      y_locf = zoo::na.locf(y, na.rm = FALSE)
    ) %>%
    ungroup()
  
  #Baseline is observed for all selected subjects, so this should remove no rows
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
  #last_obs_y is included because missingness may depend on most recently previously observed tumour burden value
  
  obs_model <- glm(
    R ~ last_obs_y + z + time + BMMTR1 + LIVERMET + B_ECOG,
    data = dat_w %>% filter(visit > 0),
    family = binomial()
  )
  
  #Estimate observation probabilities
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
    wgee_data = dat_w
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
    drgee_data = dat_dr
  )
}

#############################
#Run all methods
#############################

available_result <- fit_available_gee(prime_50)

locf_result <- fit_locf_real(prime_50)

wgee_result <- fit_wgee_real(
  prime_50,
  min_pi = 0.05
)

mice_result <- fit_mice_real(
  prime_50,
  m = 20,
  maxit = 10,
  seed = 2001
)

da_result <- fit_da_real(
  prime_50,
  n_iter = 5000,
  n_burn = 2000,
  n_thin = 2,
  n_chains = 3,
  seed = 3001
)

drgee_result <- fit_drgee_real(
  prime_50,
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

###########################
#DIAGNOSTICS
##########################


#WGEE diagnostics
#fitted observation probabilities (min is set to 0.05 by construction, mean=0.3825)
summary(wgee_result$wgee_data$pi_hat)

#inverse probability weights (mean=2.141)
summary(wgee_result$wgee_data$wgee_weight)

#how many probabilities were truncated at 0.05 (198)
sum(wgee_result$wgee_data$pi_hat == 0.05, na.rm = TRUE)



#DR-GEE diagnostics
#Check pseudo-outcomes
summary(drgee_result$drgee_data$y_dr)

#Check whether any pseudo-outcomes are extreme (-1.118,8.053)
range(drgee_result$drgee_data$y_dr, na.rm = TRUE)


#MICE diagnostics

#when ran MICE says "Number of logged events: 885", mice can still produce good results but
#logged events indicate that some of the imputation models required automatic adjustment


#DA convergence diagnostics

da_result$rhat #good values, max is <1.04

#visual check of MCMC chains
plot(da_result$da_samples) #this looks satisfactory


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
#FIGURES
#################


#Summarise missingness at each scheduled visit
visit_missingness <- prime_50 %>%
  group_by(visit, VISIT, time) %>%
  summarise(
    n_total = n(),
    n_observed = sum(R == 1, na.rm = TRUE),
    n_missing = sum(R == 0, na.rm = TRUE),
    prop_missing = mean(R == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(visit)


#missingness proportion by visit bar chart
missingness_plot <- ggplot(
  visit_missingness,
  aes(x = factor(VISIT, levels = VISIT), y = prop_missing)
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
  "missingness_prop_by_visitREAL.png",
  plot = missingness_plot,
  width = 8,
  height = 5,
  dpi = 800
)


#############################
#Lasagna plot for missingness pattern
#############################

#make sure the data are ordered correctly
lasagna_dat <- prime_50 %>%
  arrange(id, visit) %>%
  group_by(id) %>%
  mutate(
    #last observed scheduled visit for each individual
    last_observed_visit = ifelse(
      any(R == 1),
      max(visit[R == 1]),
      NA_integer_
    ),
    
    #Number of observed outcomes for each individual
    n_observed = sum(R == 1),
    
    #Indicator of intermittent missingness:
    #TRUE if a missing response is followed by a later observed response
    intermittent = any(R == 0 & rev(cummax(rev(R == 1))))
  ) %>%
  ungroup() %>%
  mutate(
    #keep visits in the correct scheduled order on the x-axis
    VISIT_ordered = factor(
      VISIT,
      levels = unique(VISIT[order(visit)])
    )
  )

#Order individuals by last observed visit and then by number of observed outcomes
#This makes patterns easier to see
id_order <- lasagna_dat %>%
  distinct(id, last_observed_visit, n_observed, intermittent) %>%
  arrange(last_observed_visit, n_observed, id) %>%
  pull(id)

lasagna_dat <- lasagna_dat %>%
  mutate(
    id_ordered = factor(id, levels = id_order)
  )

#Create lasagna plot
lasagna_plot <- ggplot(
  lasagna_dat,
  aes(x = VISIT_ordered, y = id_ordered, fill = y)
) +
  geom_tile(colour = "grey90", linewidth = 0.10) +
  
  #Observed values are shown on a cyan scale
  #Missing values are shown via gaps/white spaces
  scale_fill_gradient(
    low = "cyan1",
    high = "darkcyan",
    na.value = "white",
    name = expression(log(SLD + 1))
  ) +
  
  labs(
    x = "Scheduled visit",
    y = "Individual"
  ) +
  
  #change the theme/background to make this easier to see
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

lasagna_plot

ggsave("lasagna_plot_missingnessREAL.png",
       plot= lasagna_plot,
       width = 9,
       height = 7,
       dpi=700
       )
