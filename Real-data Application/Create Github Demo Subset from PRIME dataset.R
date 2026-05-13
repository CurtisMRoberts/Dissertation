#############################
#Create privacy-preserved GitHub subset
#############################

library(dplyr)
library(readr)
library(zoo)

set.seed(3737)

#Select 50 patients from the constructed analysis dataset
selected_ids <- prime_full %>%
  distinct(id) %>%
  slice_sample(n = 50) %>%
  pull(id)

#Create anonymous IDs
id_key <- tibble(
  id = selected_ids,
  anon_id = seq_along(selected_ids)
)

prime_full_demo_50 <- prime_full %>%
  filter(id %in% selected_ids) %>%
  left_join(id_key, by = "id") %>%
  mutate(
    #Replace original ID with anonymous ID
    id = anon_id,
    
    #add noise to the observed SLD values only
    SLD = ifelse(
      R == 1,
      pmax(
        0,
        round(SLD * exp(rnorm(n(), mean = 0, sd = 0.15)), 1)
      ),
      NA_real_
    ),
    
    #recompute transformed response from altered SLD
    y = ifelse(R == 1, log(SLD + 1), NA_real_)
  ) %>%
  select(-anon_id) %>%
  group_by(id) %>%
  arrange(visit, .by_group = TRUE) %>%
  mutate(
    baseline_y = y[visit == 0][1],
    last_obs_y = lag(zoo::na.locf(y, na.rm = FALSE)),
    last_obs_y = ifelse(visit == 0, baseline_y, last_obs_y)
  ) %>%
  ungroup() %>%
  arrange(id, visit)

#save
write_csv(prime_full_demo_50, "prime_full_demo_50.csv")
