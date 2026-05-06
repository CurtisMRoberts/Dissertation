#############
# Chapter 5 Simulation Study
#############

library(MASS)
library(dplyr)
library(tidyr)
library(geepack)
library(zoo)
library(mice)
library(rjags)
library(tidyverse)
library(VIM)
library(ggplot2)
library(scales)
library(coda)



#######################
# Function to simulate complete longitudinal Gaussian data
####################

simulatecomp_data <- function(n = 300,
                              T = 5,
                              beta = c(2, 0.4, -0.7, 0.3, -0.15),
                              rho = 0.5,
                              sigma2 = 1) {
  #visit times, T = 5 so this gives 6 repeated measurements (including an always observed baseline)
  visit <- 0:T
  n_visits <- length(visit)
  
  #exchangeable covariance matrix 
  #assumes constant variance and common correlation between any two visits
  Sigma <- matrix(rho * sigma2, nrow = n_visits, ncol = n_visits)
  diag(Sigma) <- sigma2
  
  #baseline treatment indicator
  z <- rbinom(n, size = 1, prob = 0.5)
  
  #continuous baseline covariate
  #i.e baseline disease severity 
  x <- rnorm(n, mean = 0, sd = 1)
  
  all_subjects <- list()
  
  for (i in 1:n) {
    
    #subject-specific design matrix
    X_i <- cbind(
      intercept = 1,
      time = visit,
      z = z[i],
      x = x[i],
      time_z = visit * z[i]
    )
    
    mu_i <- as.numeric(X_i %*% beta)
    
    #repeated outcomes from MVN dist
    y_i <- as.numeric(MASS::mvrnorm(n = 1, mu = mu_i, Sigma = Sigma))
    
    all_subjects[[i]] <- data.frame(
      id = i,
      visit = visit,
      time = visit,
      z = z[i],
      x = x[i],
      time_z = visit * z[i],
      y_full = y_i,
      y0 = y_i[visit == 0]
    )
  }
  
  bind_rows(all_subjects)
}

############
#MAR intermittent missingness
############

mar_missingness <- function(data_comp, missing_level = "moderate") {
  
  #these parameters control the observation probability
  #=1 means observed.
  #is MAR because depends only on observed history
  #intermittent because a subject can be missing at one visit and observed again later
  
  alpha0 <-  1.2
  alpha1 <- -0.45  #effect of last observed outcome
  alpha2 <-  0.25  #effect of treatment group z
  alpha3 <- -0.10  #effect of time t
  alpha4 <-  0.20  #effect of baseline covariate x
  
  data_out <- data_comp %>%
    group_by(id) %>%
    arrange(visit, .by_group = TRUE) %>%
    mutate(
      pi_obs = NA_real_,
      R = NA_integer_,
      y = NA_real_,
      last_obs_y = NA_real_
    ) %>%
    ungroup()
  
  split_data <- split(data_out, data_out$id)
  
  split_data <- lapply(split_data, function(d) {
    
    #baseline assumed always observed
    last_y <- d$y_full[d$visit == 0]
    
    for (r in seq_len(nrow(d))) {
      
      if (d$visit[r] == 0) {
        
        d$pi_obs[r] <- 1
        d$R[r] <- 1
        d$y[r] <- d$y_full[r]
        d$last_obs_y[r] <- last_y
        
      } else {
        
        #logistic missingness model
        #logitP(R = 1)= alpha0 + alpha1*last_obs_y + alpha2*z + alpha3*time + alpha4*x
        eta <- alpha0 +
          alpha1 * last_y +
          alpha2 * d$z[r] +
          alpha3 * d$time[r] +
          alpha4 * d$x[r]
        
        p <- 1 / (1 + exp(-eta))
        
        d$pi_obs[r] <- p
        d$R[r] <- rbinom(1, size = 1, prob = p)
        d$y[r] <- ifelse(d$R[r] == 1, d$y_full[r], NA_real_)
        d$last_obs_y[r] <- last_y
        
        #update the last observed outcome only if the current outcome is observed
        #preserves the MAR structure based on observed history
        if (d$R[r] == 1) {
          last_y <- d$y_full[r]
        }
      }
    }
    
    d
  })
  
  bind_rows(split_data) %>%
    arrange(id, visit)
}

################
#Test run for data-generating process and missingness
################

true_beta <- c(
  "(Intercept)" = 2,
  "time" = 0.4,
  "z" = -0.7,
  "x" = 0.3,
  "time:z" = -0.15
)

data_comp <- simulatecomp_data(
  n = 300,
  T = 5,
  beta = true_beta,
  rho = 0.5,
  sigma2 = 1
)

data_mar <- mar_missingness(
  data_comp,
  missing_level = "moderate"
)

#check the missingness proportion after baseline
mean(is.na(data_mar$y[data_mar$visit > 0]))


#check that the simulated missingness responds to observed history
data_mar %>%
  filter(visit > 0) %>%
  summarise(
    cor_last_y_missing = cor(last_obs_y, 1 - R, use = "complete.obs"),
    cor_x_missing = cor(x, 1 - R, use = "complete.obs"),
    cor_time_missing = cor(time, 1 - R, use = "complete.obs")
  )

#########
#LOCF
#########

fit_locf <- function(dat) {
  
  #carry last observed outcome forward for each individual
  #creates a completed outcome y_locf
  dat_locf <- dat %>%
    group_by(id) %>%
    arrange(visit, .by_group = TRUE) %>%
    mutate(
      y_locf = zoo::na.locf(y, na.rm = FALSE)
    ) %>%
    ungroup()
  
  #fit the marginal longitudinal model after imputation
  fit <- geepack::geeglm(
    y_locf ~ time * z + x,
    id = id,
    data = dat_locf,
    corstr = "exchangeable"
  )
  
  list(
    beta = coef(fit),
    se = coef(summary(fit))[, "Std.err"],
    locf_data = dat_locf,
    locf_summary = summary(fit)
  )
}

################
#WGEE
################

fit_wgee <- function(dat) {
  
  dat_w <- dat %>%
    group_by(id) %>%
    arrange(visit, .by_group = TRUE) %>%
    ungroup()
  
  #estimate observation probabilities after baseline
  #This is the missingness model used to construct inverse probability weights
  obs_model <- glm(
    R ~ last_obs_y + z + x + time,
    data = dat_w %>% filter(visit > 0),
    family = binomial()
  )
  
  dat_w <- dat_w %>%
    mutate(
      pi_hat = ifelse(
        visit == 0,
        1,
        predict(obs_model, newdata = ., type = "response")
      ),
      
      #missing rows are not used in the WGEE fit, so their weights are NA
      wgee_weight = ifelse(R == 1, 1 / pi_hat, NA_real_)
    )
  
  #fit weighted GEE using observed outcomes only
  #this is the actual outcome model of interest
  fit <- geepack::geeglm(
    y ~ time * z + x,
    id = id,
    data = dat_w %>% filter(R == 1),
    weights = wgee_weight,
    corstr = "exchangeable"
  )
  
  list(
    beta = coef(fit),
    se = coef(summary(fit))[, "Std.err"],
    obs_model = obs_model,
    wgee_data = dat_w,
    wgee_summary = summary(fit)
  )
}

#############
#MICE implementation 
#############

fit_mice <- function(dat, m = 20, maxit = 10, seed=NULL) {
  
  #convert from long to wide format (better/easier for longitudinal)
  #each individual has one row, and repeated outcomes become y0,y1,...,y5
  dat_wide <- dat %>%
    select(id, visit, z, x, y) %>%
    mutate(y_name = paste0("y", visit)) %>%
    select(id, z, x, y_name, y) %>%
    pivot_wider(
      names_from = y_name,
      values_from = y
    ) %>%
    arrange(id)
  
  #z is a treatment indicator, so treat it as categorical for imputation
  dat_wide <- dat_wide %>%
    mutate(z = factor(z))
  
  #MICE predictor matrix
  pred <- mice::make.predictorMatrix(dat_wide)
  
  #start with no predictors, then specify later what is used
  pred[,] <- 0
  
  #outcome columns to be imputed.
  y_vars <- grep("^y[0-9]+$", names(dat_wide), value = TRUE)
  
  #baseline is always observed in this simulation
  impute_vars <- setdiff(y_vars, "y0")
  
  #For each post-baseline outcome, use:
  #treatment z
  #baseline covariate x
  #other repeated outcomes
  #This reflects the chained equations specified in diss
  for (v in impute_vars) {
    pred[v, c("z", "x", setdiff(y_vars, v))] <- 1
  }
  
  #dont use id as a predictor as id is an identifier not covariate
  pred[, "id"] <- 0
  
  #specify imputation methods
  meth <- mice::make.method(dat_wide)
  meth[] <- ""
  #gaussian imputation is appropriate because the simulated outcomes are normal
  meth[impute_vars] <- "norm"
  
  #run MICE
  imp <- mice::mice(
    dat_wide,
    m = m,
    maxit = maxit,
    method = meth,
    predictorMatrix = pred,
    seed = seed,
    printFlag = FALSE
  )
  
  #fit GEE separately in each imputed dataset
  beta_list <- list()
  vcov_list <- list()
  
  for (k in 1:m) {
    
    completed_wide <- mice::complete(imp, action = k)
    
    #convert completed wide data back to long format for GEE
    completed_long <- completed_wide %>%
      pivot_longer(
        cols = all_of(y_vars),
        names_to = "visit",
        values_to = "y"
      ) %>%
      mutate(
        visit = as.numeric(gsub("y", "", visit)),
        time = visit,
        z = as.numeric(as.character(z)),
        time_z = time * z
      ) %>%
      arrange(id, visit)
    
    fit_k <- geepack::geeglm(
      y ~ time * z + x,
      id = id,
      data = completed_long,
      corstr = "exchangeable"
    )
    
    beta_list[[k]] <- coef(fit_k)
    vcov_list[[k]] <- fit_k$geese$vbeta
  }
  
  #Rubin's rules for pooling
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


#################
#Data Augmentation 
#################

fit_da <- function(dat,
                   n_iter = 5000,
                   n_burn = 2000,
                   n_thin = 2,
                   n_chains = 3,
                   seed = NULL) {
  
  #missing y values are left as NA. JAGS then treats these missing outcomes as unknown latent
  #quantities and samples them inside the MCMC algorithm
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #convert to wide format so that each subject has one row and the repeated outcomes are stored y0,y1,...
  dat_wide <- dat %>%
    select(id, visit, z, x, y) %>%
    mutate(y_name = paste0("y", visit)) %>%
    select(id, z, x, y_name, y) %>%
    pivot_wider(
      names_from = y_name,
      values_from = y
    ) %>%
    arrange(id)
  
  #identify the repeated outcome columns in the correct order
  #should give c("y0", "y1", ..., "y5")
  y_vars <- paste0("y", 0:5)
  
  #Response matrix for JAGS.
  #rows correspond to subjects and columns correspond to visits
  Y <- as.matrix(dat_wide[, y_vars])
  
  #number of individuals
  n_id <- nrow(Y)
  
  
  #Subject-level covariates generated at baseline and remain fixed over visits
  z <- as.numeric(dat_wide$z)
  x <- as.numeric(dat_wide$x)
  
  #visit times
  time <- 0:5
  
  #JAGS needs a precision matrix, not a covariance matrix
  #construct the covariance matrix in JAGS, then invert it
  jags_data <- list(
    n_id = n_id,
    Y = Y,
    z = z,
    x = x,
    time = time
  )
  
  model_string <- "
  model {

    ####################
    #Likelihood
    ####################

    for (i in 1:n_id) {

      # Marginal mean model.
      # time[1] = 0, time[2] = 1, ..., time[6] = 5.
      for (j in 1:6) {
        mu[i, j] <- beta[1] +
                    beta[2] * time[j] +
                    beta[3] * z[i] +
                    beta[4] * x[i] +
                    beta[5] * time[j] * z[i]
      }

      #baseline
      Y[i, 1] ~ dnorm(mu[i, 1], tau1)
      resid[i, 1] <- Y[i, 1] - mu[i, 1]

      #Visit 1
      cond_mean[i, 2] <- mu[i, 2] +
                         c2 * resid[i, 1]

      Y[i, 2] ~ dnorm(cond_mean[i, 2], tau2)
      resid[i, 2] <- Y[i, 2] - mu[i, 2]

      #Visit 2
      cond_mean[i, 3] <- mu[i, 3] +
                         c3 * (resid[i, 1] + resid[i, 2])

      Y[i, 3] ~ dnorm(cond_mean[i, 3], tau3)
      resid[i, 3] <- Y[i, 3] - mu[i, 3]

      #Visit 3
      cond_mean[i, 4] <- mu[i, 4] +
                         c4 * (resid[i, 1] + resid[i, 2] + resid[i, 3])

      Y[i, 4] ~ dnorm(cond_mean[i, 4], tau4)
      resid[i, 4] <- Y[i, 4] - mu[i, 4]

      #Visit 4
      cond_mean[i, 5] <- mu[i, 5] +
                         c5 * (resid[i, 1] + resid[i, 2] + resid[i, 3] + resid[i, 4])

      Y[i, 5] ~ dnorm(cond_mean[i, 5], tau5)
      resid[i, 5] <- Y[i, 5] - mu[i, 5]

      #Visit 5
      cond_mean[i, 6] <- mu[i, 6] +
                         c6 * (resid[i, 1] + resid[i, 2] + resid[i, 3] +
                               resid[i, 4] + resid[i, 5])

      Y[i, 6] ~ dnorm(cond_mean[i, 6], tau6)
      resid[i, 6] <- Y[i, 6] - mu[i, 6]
    }


    #############################
    #Conditional MVN coefficients
    ###########################

    #These conditional distributions are equivalent to the 6x6 exchangable MVN model
    #but avoids using dmnorm() on a partially observed vector (as that causes an error in jags)
    
    #baseline variance
    tau1 <- 1 / pow(sigma, 2)

    #conditional coefficient and precision for visit 1 given baseline
    c2 <- rho
    var2 <- pow(sigma, 2) * (1 - pow(rho, 2))
    tau2 <- 1 / var2

    #conditional coefficient and precision for visit 2 given previous outcomes
    c3 <- rho / (1 + rho)
    var3 <- pow(sigma, 2) * (1 - (2 * pow(rho, 2)) / (1 + rho))
    tau3 <- 1 / var3

    #conditional coefficient and precision for visit 3 given previous outcomes
    c4 <- rho / (1 + 2 * rho)
    var4 <- pow(sigma, 2) * (1 - (3 * pow(rho, 2)) / (1 + 2 * rho))
    tau4 <- 1 / var4

    #conditional coefficient and precision for visit 4 given previous outcomes
    c5 <- rho / (1 + 3 * rho)
    var5 <- pow(sigma, 2) * (1 - (4 * pow(rho, 2)) / (1 + 3 * rho))
    tau5 <- 1 / var5

    #conditional coefficient and precision for visit 5 given previous outcomes
    c6 <- rho / (1 + 4 * rho)
    var6 <- pow(sigma, 2) * (1 - (5 * pow(rho, 2)) / (1 + 4 * rho))
    tau6 <- 1 / var6


    #############
    #Priors
    ############

    #normal priors for the regression coefficients
    #a precision of 0.01 corresponds to a variance of 100
    for (p in 1:5) {
      beta[p] ~ dnorm(0, 0.01)
    }

    #weakly informative prior for the sd
    sigma ~ dunif(0, 10)

    #prior for correlation parameter, set slightly below bounds for valid covariance matrix 
    rho ~ dunif(-0.199, 0.999)
  }
  "
  
  #Chain-specific initial values.
  #these are chosen to avoid identical initial values across chains
  #the initial values are also kept inside the parameter space
  inits <- function(chain = 1) {
    list(
      beta = rnorm(5, mean = 0, sd = 0.2),
      sigma = runif(1, min = 0.7, max = 1.5),
      rho = runif(1, min = 0.1, max = 0.8),
      .RNG.name = "base::Mersenne-Twister", #chosen random number generator
      .RNG.seed = ifelse(is.null(seed), sample.int(1e6, 1), seed + chain) #makes results reproducible 
    )
  }
  
  init_list <- lapply(seq_len(n_chains), inits)
  
  #Compile the model
  da_model <- rjags::jags.model(
    file = textConnection(model_string),
    data = jags_data,
    inits = init_list,
    n.chains = n_chains,
    n.adapt = 1000,
    quiet = TRUE
  )
  
  #the burn-in period
  update(da_model, n.iter = n_burn, progress.bar = "none")
  
  #Drawing posterior samples
  da_samples <- rjags::coda.samples(
    model = da_model,
    variable.names = c("beta", "sigma", "rho"),
    n.iter = n_iter,
    thin = n_thin,
    progress.bar = "none"
  )
  
  #Convert MCMC output to a matrix for posterior summaries
  sample_mat <- as.matrix(da_samples)
  
  #Extract beta samples 
  beta_names_jags <- paste0("beta[", 1:5, "]")
  beta_samples <- sample_mat[, beta_names_jags, drop = FALSE]
  
  #Posterior means are used as point estimates
  beta_hat <- colMeans(beta_samples)
  
  #Posterior standard deviations are used as uncertainty estimates
  se_hat <- apply(beta_samples, 2, sd)
  
  #do this to rename regression parameters to match the other methods
  term_names <- c("(Intercept)", "time", "z", "x", "time:z")
  names(beta_hat) <- term_names
  names(se_hat) <- term_names
  
  #Convergence diagnostics
  #returns Rhat for beta, sigma, and rho to help quantify convergence
  rhat <- tryCatch(
    {
      diag <- coda::gelman.diag(da_samples, autoburnin = FALSE, multivariate = FALSE)
      rhat_out <- diag$psrf[, "Point est."]
      rhat_out
    },
    error = function(e) {
      rep(NA_real_, length(c(beta_names_jags, "rho", "sigma")))
    }
  )
  
  #need to rename beta Rhat values so they match the regression terms
  if (!is.null(rhat)) {
    names(rhat) <- gsub("beta\\[1\\]", "(Intercept)", names(rhat))
    names(rhat) <- gsub("beta\\[2\\]", "time", names(rhat))
    names(rhat) <- gsub("beta\\[3\\]", "z", names(rhat))
    names(rhat) <- gsub("beta\\[4\\]", "x", names(rhat))
    names(rhat) <- gsub("beta\\[5\\]", "time:z", names(rhat))
  }
  
  list(
    beta = beta_hat,
    se = se_hat,
    da_samples = da_samples,
    beta_samples = beta_samples,
    rhat = rhat,
    model = da_model,
    da_data_wide = dat_wide,
    Y_matrix = Y
  )
}


################
# DR-GEE 
################

fit_drgee <- function(dat,
                      min_pi = 0.05,
                      corstr = "exchangeable") {

  
  
  #the augmented pseudo-outcome form is specified in dissertation - essentially creats Y that acts like full-data outcome
  #if the missingness model is correct, the inverse probability part corrects bias
  #if the outcome model is correct, the augmentation part corrects bias
  #hence it is doubly robust
  
  dat_dr <- dat %>%
    group_by(id) %>%
    arrange(visit, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(
      z = as.numeric(z),
      time_z = time * z
    )
  
  
  #missingness model
  #this is the same observation model used in WGEE
  pi_model <- glm(
    R ~ last_obs_y + z + x + time,
    data = dat_dr %>% filter(visit > 0),
    family = binomial()
  )
  
  dat_dr <- dat_dr %>%
    mutate(
      pi_hat_raw = ifelse(
        visit == 0,
        1,
        predict(pi_model, newdata = ., type = "response")
      ),
      
      #if pi_hat is extremely small, the inverse weight becomes huge
      #this can make the estimating equation unstable so truncating 
      #very small probabilities for numerical stability is done
      pi_hat = ifelse(
        visit == 0,
        1,
        pmax(pi_hat_raw, min_pi)
      ),
      
      ipw_weight = ifelse(R == 1, 1 / pi_hat, NA_real_)
    )
  
  
  #fit the outcome model for augmentation
  #the augmentation model estimates the expected outcome at each visit,
  #conditional on observed information
  
  outcome_model <- lm(
    y ~ time * z + x + last_obs_y,
    data = dat_dr %>% filter(R == 1, visit > 0)
  )
  
  dat_dr <- dat_dr %>%
    mutate(
      m_hat = ifelse(
        visit == 0,
        y,
        predict(outcome_model, newdata = .)
      )
    )
  
  
  #construct the DR pseudo-outcome
  #for missing rows, y is NA, so we replace y by 0
  #this works because R = 0 for missing rows, so the observed outcome term vanishes
  
  dat_dr <- dat_dr %>%
    mutate(
      y_obs_for_formula = ifelse(R == 1, y, 0),
      
      y_dr = (R / pi_hat) * y_obs_for_formula +
        (1 - R / pi_hat) * m_hat
    )
  

  #fit GEE to the DR pseudo-outcome
  #now fit the same marginal mean model as the other methods
  #but the response is now y_dr rather than y
  
  fit <- geepack::geeglm(
    y_dr ~ time * z + x,
    id = id,
    data = dat_dr,
    corstr = corstr
  )
  
  ################
  #return estimates
  ################
  
  beta_hat <- coef(fit)
  se_hat <- coef(summary(fit))[, "Std.err"]
  
  list(
    beta = beta_hat,
    se = se_hat,
    drgee_summary = summary(fit),
    drgee_data = dat_dr,
    pi_model = pi_model,
    outcome_model = outcome_model,
    pi_summary = summary(pi_model),
    outcome_summary = summary(outcome_model)
  )
}

#######################
#Simulation loop for replications 
######################
true_beta <- c(
  "(Intercept)" = 2,
  "time" = 0.4,
  "z" = -0.7,
  "x" = 0.3,
  "time:z" = -0.15
)

results_list <- list()
#store DA convergence diagnostic for each simulation, this is a diagnostic not an estimate 
da_rhat_list <-list()
for (s in 1:200) {
  
  set.seed(1234+ s)
  
  data_comp <- simulatecomp_data(
    n = 300,
    T = 5,
    beta = true_beta,
    rho = 0.5,
    sigma2 = 1
  )
  
  data_mar <- mar_missingness(data_comp)
  
  fit_complete <- geepack::geeglm(
    y_full ~ time * z + x,
    id = id,
    data = data_comp,
    corstr = "exchangeable"
  )
  
  locf_result <- fit_locf(data_mar)
  wgee_result <- fit_wgee(data_mar)
  mice_result <- fit_mice(data_mar, m = 20, maxit = 10, seed = 1000 + s)
  
  da_result <- fit_da(
    data_mar,
    n_iter = 5000,
    n_burn = 2000,
    n_thin = 2,
    n_chains = 3,
    seed = 2000+s
  )
  #store Rhat for parameters 
  da_rhat_list[[s]] <- data.frame(
    sim = s,
    parameter = names(da_result$rhat),
    rhat = as.numeric(da_result$rhat)
  )
  
  drgee_result <- fit_drgee(data_mar)
  
  results_list[[s]] <- bind_rows(
    data.frame(sim = s, method = "Complete", term = names(coef(fit_complete)),
               estimate = as.numeric(coef(fit_complete)),
               se = as.numeric(coef(summary(fit_complete))[, "Std.err"])),
    
    data.frame(sim = s, method = "LOCF", term = names(locf_result$beta),
               estimate = as.numeric(locf_result$beta),
               se = as.numeric(locf_result$se)),
    
    data.frame(sim = s, method = "WGEE", term = names(wgee_result$beta),
               estimate = as.numeric(wgee_result$beta),
               se = as.numeric(wgee_result$se)),
    
    data.frame(sim = s, method = "MICE", term = names(mice_result$beta),
               estimate = as.numeric(mice_result$beta),
               se = as.numeric(mice_result$se)),
    
    data.frame(sim = s, method = "DA", term = names(da_result$beta),
               estimate = as.numeric(da_result$beta),
               se = as.numeric(da_result$se)),
    
    data.frame(sim = s, method = "DRGEE", term = names(drgee_result$beta),
               estimate = as.numeric(drgee_result$beta),
               se = as.numeric(drgee_result$se))
  )
  
  cat("Completed simulation", s, "of", 200, "\n")
}

sim_results <- bind_rows(results_list)


true_values <- data.frame(
  term = names(true_beta),
  true = as.numeric(true_beta)
)

sim_summary <- sim_results %>%
  left_join(true_values, by = "term") %>%
  group_by(method, term) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - true, na.rm = TRUE),
    empirical_sd = sd(estimate, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    rmse = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
    coverage = mean(
      estimate - 1.96 * se <= true &
        estimate + 1.96 * se >= true,
      na.rm = TRUE
    ),
    n_success = sum(!is.na(estimate)),
    .groups = "drop"
  )

sim_summary

#combine all DA rhat values into one
da_rhat_results <- bind_rows(da_rhat_list)

#summarise convergence diagnostics across simulation replications
da_rhat_summary <- da_rhat_results %>%
  group_by(parameter) %>%
  summarise(
    mean_rhat = mean(rhat, na.rm = TRUE),
    median_rhat = median(rhat, na.rm = TRUE),
    max_rhat = max(rhat, na.rm = TRUE),
    n_available = sum(!is.na(rhat)),
    .groups = "drop"
  )

da_rhat_summary


###########
#FIGURES
##########

#missingness proportion at each visit 
missing_by_visit <- data_mar %>%
  group_by(visit) %>%
  summarise(
    n = n(),
    n_missing = sum(is.na(y)),
    prop_missing = mean(is.na(y)),
    .groups = "drop"
  )


p_missing_visit <- ggplot(missing_by_visit, aes(x = factor(visit), y = prop_missing)) +
  geom_col(width = 0.65, fill = "grey35") +
  geom_text(
    aes(label = percent(prop_missing, accuracy = 1)),
    vjust = -0.35,
    size = 3.8
  ) +
  scale_y_continuous(
    limits = c(0, 0.75),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = "Visit",
    y = "Proportion missing"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

p_missing_visit

ggsave(
  "C:/Users/curti/Downloads/missingness_by_visit.png",
  width = 8,
  height = 6,
  dpi = 300
)

#missingness heatmap

missing_heatmap <- data_mar %>%
  filter(id <= 75) %>%
  mutate(
    status = ifelse(is.na(y), "Missing", "Observed"),
    id = factor(id),
    visit = factor(visit)
  )

ggplot(missing_heatmap, aes(x = visit, y = id, fill = status)) +
  geom_tile(colour = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Observed" = "white", "Missing" = "darkcyan")
  ) +
  labs(
    x="Visit",
    y="Subject",
    fill="Status"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

ggsave(
  "C:/Users/curti/Downloads/missingness_heatmap.png",
  width = 8,
  height = 6,
  dpi = 300
)

#coverage plot

#coverage results from the updated simulation summary
coverage_results <- data.frame(
  method = c(
    rep("Complete", 5),
    rep("DA", 5),
    rep("DR-GEE", 5),
    rep("LOCF", 5),
    rep("MICE", 5),
    rep("WGEE", 5)
  ),
  parameter = rep(
    c("beta_0", "beta_time", "beta_time:z", "beta_x", "beta_z"),
    times = 6
  ),
  coverage = c(
    #Complete
    0.970, 0.965, 0.970, 0.940, 0.965,
    #DA
    0.960, 0.970, 0.985, 0.945, 0.965,
    #DR-GEE
    0.945, 0.930, 0.960, 0.950, 0.965,
    #LOCF
    0.965, 0.000, 0.195, 0.945, 0.955,
    #MICE
    0.950, 0.935, 0.950, 0.945, 0.955,
    #WGEE
    0.905, 0.925, 0.960, 0.935, 0.960
  )
)

#order methods and parameters
coverage_results$method <- factor(
  coverage_results$method,
  levels = c("Complete", "MICE", "DR-GEE", "WGEE", "DA", "LOCF")
)

coverage_results$parameter <- factor(
  coverage_results$parameter,
  levels = c("beta_0", "beta_time", "beta_time:z", "beta_x", "beta_z")
)

coverage_plot <- ggplot(
  coverage_results,
  aes(x = parameter, y = coverage, group = method)
) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ method, ncol = 3) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 0.95, 1.00)
  ) +
  scale_x_discrete(
    labels = c(
      "beta_0" = expression(beta[0]),
      "beta_time" = expression(beta[time]),
      "beta_time:z" = expression(beta[time:z]),
      "beta_x" = expression(beta[x]),
      "beta_z" = expression(beta[z])
    )
  ) +
  labs(
    x="Parameter",
    y="Empirical coverage"
  ) +
  theme_bw()

coverage_plot

ggsave(
  "C:/Users/curti/Downloads/coverage_plot.png",
  width = 8,
  height = 6,
  dpi = 500
)
