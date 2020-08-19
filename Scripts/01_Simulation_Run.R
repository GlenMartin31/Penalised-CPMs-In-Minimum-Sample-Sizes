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
## This script runs the simulations across all scenarios: runs in parallel using furrr
####-----------------------------------------------------------------------------------------

#Load the simulation functions
source(here::here("Scripts", "Simulation_Functions.R"))
library(tidyverse)
library(furrr)

# Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  P = 10,
  RhoX = c(0, 0.5),
  Y_prev = c(0.2, 0.5),
  R2_based_on_maxR2 = c(TRUE, FALSE)
)

# number of repeats per scenario
n_rep <- 100

#Set up parallel process
plan(multiprocess, workers = (availableCores() - 1))
sims_allXpredict <- sims_parameters %>%
  mutate(results = future_pmap(list(n_iter = n_rep, 
                                    P = P, 
                                    RhoX = RhoX, 
                                    Y_prev = Y_prev, 
                                    R2_based_on_maxR2 = R2_based_on_maxR2),
                               simulation_nruns_fnc,
                               beta_true = c(rep(log(1.5), 6), rep(log(2), 2), rep(log(3), 2)),
                               .progress = TRUE,
                               .options = future_options(seed = as.integer(698137)))
         )
write_rds(sims_allXpredict, path = here::here("Data", "sims_allXpredict.RDS"))


#Set up parallel process
plan(multiprocess, workers = (availableCores() - 1))
sims_halfXpredict <- sims_parameters %>%
  mutate(results = future_pmap(list(n_iter = n_rep, 
                                    P = P, 
                                    RhoX = RhoX, 
                                    Y_prev = Y_prev, 
                                    R2_based_on_maxR2 = R2_based_on_maxR2),
                               simulation_nruns_fnc,
                               beta_true = c(rep(log(1.5), 3), log(2), log(3), rep(0, 5)),
                               .progress = TRUE,
                               .options = future_options(seed = as.integer(91605)))
         )
write_rds(sims_halfXpredict, path = here::here("Data", "sims_halfXpredict.RDS"))
