# #######################################################################################################################

# Author of code: Glen P. Martin.

# This is code for a simulation study presented in a manuscript entitled: 
# Developing Clinical Prediction Models Using Data that Adheres to
# Minimum Sample Size Criteria: the importance of penalization methods and quantifying bootstrap variability
# Authors:
#   Glen P. Martin
#   Richard Riley
#   Gary S. Collins
#   Matthew Sperrin


# #######################################################################################################################

####-----------------------------------------------------------------------------------------
## This script summarises the simulation results, and manipulates the data in a format ready
## for plotting.
####-----------------------------------------------------------------------------------------
library(tidyverse)

sims_allXpredict <- read_rds(here::here("Data", "sims_allXpredict.RDS"))
sims_halfXpredict <- read_rds(here::here("Data", "sims_halfXpredict.RDS"))

###-------------------------------------------------------------------------------------
## Combine Simulation scenarios
###-------------------------------------------------------------------------------------
sims_all <- sims_allXpredict %>%
  mutate("Beta" = "allX") %>%
  bind_rows(sims_halfXpredict %>%
              mutate("Beta" = "halfX")) %>%
  mutate("Beta" = factor(Beta, levels = c("allX", "halfX")),
         
         results = map(results, as_tibble)) %>%
  unnest(cols = c(results)) 


#Create a simulation scenario indicator:
sims_all <- sims_all %>%
  arrange(RhoX, Y_prev, R2_based_on_maxR2, Beta) %>%
  mutate("SimulationScenario" = group_indices(., P, RhoX, Y_prev, R2_based_on_maxR2, Beta)) %>%
  select(SimulationScenario, P, RhoX, Y_prev, R2_based_on_maxR2, Beta, everything()) %>%
  ungroup()

# write_rds(sims_all, path = here::here("Data", "sims_AllResults.RDS"))




###-------------------------------------------------------------------------------------
## Summarise the performance metrics that were generated
###-------------------------------------------------------------------------------------

mean_summary_sims <- sims_all %>%
  group_by(SimulationScenario) %>%
  summarise_at(.vars = vars(R2_mle:R2_ridge,
                            starts_with("CITL_"),
                            starts_with("CalSlope_"),
                            starts_with("AUC_")),
               .funs = list("mean" = mean,
                            "sd" = sd)) %>%
  ungroup() %>%
  #Turn into long format:
  pivot_longer(cols = -SimulationScenario, 
               names_to = "metric_model_quantity", 
               values_to = "value") %>%
  #Separate the metric_model_quantity variable into separate colums based on regex.:
  extract(metric_model_quantity, 
          into = c("PerformanceMetric", "Model", "Quantity"), 
          "([A-Za-z2]+)_([A-Za-z]+)_([A-Za-z]+)") %>%
  #Place the SD estimates 'next to' the mean estimates for each model, scenario and performance metric combination:
  pivot_wider(id_cols = c("SimulationScenario", "Model", "PerformanceMetric"),
              names_from = "Quantity",
              values_from = "value")


###-------------------------------------------------------------------------------------
## Summarise the distribution of the performance metrics across all simulations
###-------------------------------------------------------------------------------------
CITL_IQR <- sims_all %>%
  select(SimulationScenario,
         starts_with("CITL_")) %>% 
  pivot_longer(cols = c(starts_with("CITL_")), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
  group_by(Model, SimulationScenario) %>%
  summarise("LQ" = quantile(value, probs = 0.25),
            "UQ" = quantile(value, probs = 0.75),
            "Width" = UQ - LQ)


CalSlope_IQR <- sims_all %>%
  select(SimulationScenario,
         starts_with("CalSlope_")) %>% 
  pivot_longer(cols = c(starts_with("CalSlope_")), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
  group_by(Model, SimulationScenario) %>%
  summarise("LQ" = quantile(value, probs = 0.25),
            "UQ" = quantile(value, probs = 0.75),
            "Width" = UQ - LQ)

