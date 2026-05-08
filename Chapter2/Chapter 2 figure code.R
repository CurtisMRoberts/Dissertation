rm(list=ls())

#packages
library(ggplot2)

#Create the dataset
data <- data.frame(
  IQ = c(78, 84, 84, 85, 87, 91, 92, 94, 94, 96, 99, 105, 105, 106, 108, 112, 113, 115, 118, 134),
  JobPerformance = c(9, 13, 10, 8, 7, 7, 9, 9, 11, 7, 7, 10, 11, 15, 10, 10, 12, 14, 16, 12)
)

#Suppose JobPerformance is missing for IQ < 99
data$Observed <- ifelse(data$IQ < 99, FALSE, TRUE)


###############################
#Complette data scatterplot
###############################

ggplot(data, aes(x = IQ, y = JobPerformance)) +
  geom_point(size = 3) +
  labs(
    title = "Complete-data scatterplot of IQ and job performance",
    x = "IQ",
    y = "Job Performance",
    color = "Status"
  ) +
  theme_minimal(base_size = 14)


################################
#Listwise deletion scatterplot
################################

#Keep only observed cases (IQ >= 99)
observed_data <- subset(data, IQ >= 99)

#Plot with original IQ range
ggplot(observed_data, aes(x = IQ, y = JobPerformance)) +
  geom_point(size = 3) +
  labs(
    x = "IQ",
    y = "Job Performance"
  ) +
  scale_x_continuous(limits = range(data$IQ))


###############################
#Mean imputation scatterplot 
###############################
rm(list=ls())


#original complete dataset
data <- data.frame(
  IQ = c(78, 84, 84, 85, 87, 91, 92, 94, 94, 96,
         99, 105, 105, 106, 108, 112, 113, 115, 118, 134),
  JobPerformance = c(9, 13, 10, 8, 7, 7, 9, 9, 11, 7,
                     7, 10, 11, 15, 10, 10, 12, 14, 16, 12)
)

#job performance is missing for individuals with IQ < 99
data$Observed <- data$IQ >= 99

#Split the data into observed and missing subsets
observed <- subset(data, Observed == TRUE)
missing  <- subset(data, Observed == FALSE)

#Calc the mean of the observed job performance scores only (value used to replace missing values)
observed_mean <- mean(observed$JobPerformance)

#Create a new variable containing the imputed job performance values,
#if the value is observed, keep the original value. if missing replace it with observed mean
data$JobPerformance_Imputed <- ifelse(
  data$Observed == TRUE,
  data$JobPerformance,
  observed_mean
)

#splits the data into observed and imputed for colouring in the plot
data$Status <- ifelse(data$Observed == TRUE, "Observed", "Imputed")

#mean imputation scatterplot
mean_imp_scatter <- ggplot(data, aes(x = IQ, y = JobPerformance_Imputed, colour = Status)) +
  
  geom_point(size = 3) +
  
  #adds a dashed horizontal line showing the observed mean used for imputation
  geom_hline(yintercept = observed_mean, linetype = "dashed") +
  
  labs(
    title = "Mean Imputation",
    subtitle = "Missing values are replaced by the observed mean",
    x = "IQ",
    y = "Job Performance",
    colour = "Data Status"
  ) +
  
  #observed points are black and imputed points are red
  scale_colour_manual(values = c("Observed" = "black", "Imputed" = "#e31a1c")) +
  
  #Keep the original complete-data IQ range 
  scale_x_continuous(limits = range(data$IQ)) 

ggsave(
  filename = "mean_imputation_plot.png",
  plot = mean_imp_scatter,
  width = 6,
  height = 4,
)
  

###################################
#Stochastic regression imputation
###################################
rm(list=ls())
#Original dataset again
data <- data.frame(
  IQ = c(78, 84, 84, 85, 87, 91, 92, 94, 94, 96, 99, 105, 105, 106, 108, 112, 113, 115, 118, 134),
  JobPerformance = c(9, 13, 10, 8, 7, 7, 9, 9, 11, 7, 7, 10, 11, 15, 10, 10, 12, 14, 16, 12)
)

#Split observed and missing (simulate missing for IQ < 99)
observed <- subset(data, IQ >= 99)
missing <- subset(data, IQ < 99)

#Fit regression on observed data
fit <- lm(JobPerformance ~ IQ, data = observed)

# Stochastic regression imputation: predicted + residual noise
set.seed(123)  # for reproducibility
residual_sd <- sd(residuals(fit))  # standard deviation of residuals
missing$JobPerformance <- predict(fit, newdata = missing) + rnorm(nrow(missing), mean = 0, sd = residual_sd)

# Combine data for plotting
imputed_data <- rbind(observed, missing)
imputed_data$Status <- ifelse(imputed_data$IQ >= 99, "Observed", "Imputed")

# Plot
ggplot(imputed_data, aes(x = IQ, y = JobPerformance, color = Status)) +
  geom_point(size = 3) +
  geom_smooth(data = observed, aes(x = IQ, y = JobPerformance), method = "lm", se = FALSE, color = "black") +
  labs(
    title = "Stochastic Regression Imputation",
    subtitle = "Imputed points scatter around regression line",
    x = "IQ",
    y = "Job Performance",
    color = "Data Status"
  ) +
  scale_color_manual(values = c("Observed" = "black", "Imputed" = "#e31a1c")) +
  scale_x_continuous(limits = range(data$IQ)) 


###############################################################
# MULTIPLE IMPUTATION ILLUSTRATION
###############################################################
rm(list=ls())


#original complete dataset
data <- data.frame(
  IQ = c(78, 84, 84, 85, 87, 91, 92, 94, 94, 96,
         99, 105, 105, 106, 108, 112, 113, 115, 118, 134),
  JobPerformance = c(9, 13, 10, 8, 7, 7, 9, 9, 11, 7,
                     7, 10, 11, 15, 10, 10, 12, 14, 16, 12)
)

#Job performance is treated as missing for individuals with IQ < 99
data$Observed <- data$IQ >= 99

#Split the data into observed and missing subsets
observed <- subset(data, Observed == TRUE)
missing  <- subset(data, Observed == FALSE)

#Fit a regression model using only the observed cases
fit <- lm(JobPerformance ~ IQ, data = observed)

#Estimate the residual standard deviation from the fitted regression model
#used to add random residual variation to the imputed values
residual_sd <- sd(residuals(fit))

#seed so that the random imputations are reproducible
set.seed(123)

#the number of imputed datasets to generate
K <- 10

#Generate K completed datasets
mi_list <- lapply(1:K, function(k) {
  
  #copy the original data for the kth imputed dataset
  temp <- data
  
  #Stores the imputation number
  temp$Imputation <- paste0("Imputation ", k)
  
  #Create a new variable that will contain observed and imputed values
  temp$JobPerformance_Imputed <- temp$JobPerformance
  
  #predict job performance for the missing cases using the regression model
  predicted_values <- predict(fit, newdata = temp[temp$Observed == FALSE, ])
  
  #Add random residual noise to the predicted values
  #gives a different plausible set of imputed values for each dataset
  temp$JobPerformance_Imputed[temp$Observed == FALSE] <-
    predicted_values + rnorm(
      n = sum(temp$Observed == FALSE),
      mean = 0,
      sd = residual_sd
    )
  
  #Create a variable for colouring points in the plot
  temp$Status <- ifelse(temp$Observed == TRUE, "Observed", "Imputed")
  
  #Returns the completed dataset
  return(temp)
})

#Combine K completed datasets into one dataset for plotting
mi_data <- do.call(rbind, mi_list)

# For the observed points, keep only one copy, otherwise the observed points would be plotted K times on top of each other

observed_once <- subset(mi_data, Status == "Observed" & Imputation == "Imputation 1")

#Keep all imputed points from all K imputations
imputed_all <- subset(mi_data, Status == "Imputed")

#MI scatterplot
MI_scatter<-ggplot() +
  
  #Plot all imputed values from all K imputations
  #Alpha makes the points semi-transparent, helps show overlap
  geom_point(
    data = imputed_all,
    aes(x = IQ, y = JobPerformance_Imputed, colour = Status),
    size = 2,
    alpha = 0.45
  ) +
  
  #Plot the observed values once
  geom_point(
    data = observed_once,
    aes(x = IQ, y = JobPerformance_Imputed, colour = Status),
    size = 3
  ) +
  
  #Add the fitted regression line from the observed data
  geom_smooth(
    data = observed,
    aes(x = IQ, y = JobPerformance),
    method = "lm",
    se = FALSE,
    colour = "black"
  ) +
  
  
  labs(
    title = "Multiple Imputation Illustration",
    subtitle = "Several plausible imputed values are generated for each missing case",
    x = "IQ",
    y = "Job Performance",
    colour = "Data Status"
  ) +
  
  #observed points black and imputed points red
  scale_colour_manual(values = c("Observed" = "black", "Imputed" = "#e31a1c")) +
  
  
  scale_x_continuous(limits = range(data$IQ)) 
  

ggsave(
  filename = "MI_plot.png",
  plot = MI_scatter,
  width = 6,
  height = 4,
)
