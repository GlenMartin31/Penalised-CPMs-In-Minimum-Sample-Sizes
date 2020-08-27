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

sims_all <- read_rds(here::here("Data", "simulation_results_all.RDS")) %>%
  arrange(Simulation_Scenario) %>%
  rename("Beta" = "beta_true") %>%
  mutate("Beta" = factor(Beta, levels = c("AllX", "HalfX")),
         
         Model = fct_relevel(fct_recode(Model
                                        , "MLE" = "MLE"
                                        , "Uniform closed-form" = "Uniform"
                                        , "Uniform bootstrap" = "Bootstrap"
                                        , "Firths" = "Firth"
                                        , "LASSO" = "LASSO"
                                        , "Ridge" = "Ridge"
         ),
         "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
         "LASSO", "Ridge"))



###-------------------------------------------------------------------------------------
## Summarise the performance metrics that were generated
###-------------------------------------------------------------------------------------

mean_summary_sims <- sims_all %>%
  group_by(Simulation_Scenario, Model) %>%
  summarise_at(.vars = vars(R2,
                            starts_with("CITL"),
                            starts_with("CalSlope"),
                            starts_with("AUC")),
               .funs = list("mean" = mean,
                            "BetweenSD" = sd)) %>%
  ungroup() %>%
  #Turn into long format:
  pivot_longer(cols = c(-Simulation_Scenario, -Model), 
               names_to = "metric_quantity", 
               values_to = "value") %>%
  #Separate the metric_model_quantity variable into separate colums based on regex.:
  extract(metric_quantity, 
          into = c("PerformanceMetric", "Quantity"), 
          "([A-Za-z2]+)_([A-Za-z]+)") %>%
  #Place the SD estimates 'next to' the mean estimates for each model, scenario and performance metric combination:
  pivot_wider(id_cols = c("Simulation_Scenario", "Model", "PerformanceMetric"),
              names_from = "Quantity",
              values_from = "value")


###-------------------------------------------------------------------------------------
## Summarise the distribution of the performance metrics across all simulations
###-------------------------------------------------------------------------------------
CITL_IQR <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         CITL) %>% 
  group_by(Model, Simulation_Scenario) %>%
  summarise("LQ" = quantile(CITL, probs = 0.25),
            "UQ" = quantile(CITL, probs = 0.75),
            "Width" = UQ - LQ,
            .groups = "drop") %>%
  ungroup()


CalSlope_IQR <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         CalSlope) %>% 
  group_by(Model, Simulation_Scenario) %>%
  summarise("LQ" = quantile(CalSlope, probs = 0.25),
            "UQ" = quantile(CalSlope, probs = 0.75),
            "Width" = UQ - LQ,
            .groups = "drop") %>%
  ungroup()

