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

#Set up paralell process
plan(multiprocess)

set.seed(698137)
sims_allXpredict <- sims_parameters %>%
  mutate(results = future_pmap(list(n_iter = n_rep, 
                                    P = P, 
                                    RhoX = RhoX, 
                                    Y_prev = Y_prev, 
                                    R2_based_on_maxR2 = R2_based_on_maxR2),
                               simulation_nruns_fnc,
                               beta_true = c(rep(0.2, 6), rep(0.5, 2), rep(0.8, 2)),
                               .progress = TRUE
                               )
         )

write_rds(sims_allXpredict, path = here::here("Data", "sims_allXpredict.RDS"))

set.seed(698138)
sims_halfXpredict <- sims_parameters %>%
  mutate(results = future_pmap(list(n_iter = n_rep, 
                                    P = P, 
                                    RhoX = RhoX, 
                                    Y_prev = Y_prev, 
                                    R2_based_on_maxR2 = R2_based_on_maxR2),
                               simulation_nruns_fnc,
                               beta_true = c(0.2, 0.2, 0.2, 0.5, 0.8, rep(0, 5)),
                               .progress = TRUE
                               )
         )

write_rds(sims_halfXpredict, path = here::here("Data", "sims_halfXpredict.RDS"))
