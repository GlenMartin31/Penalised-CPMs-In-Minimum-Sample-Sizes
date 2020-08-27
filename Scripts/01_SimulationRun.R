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

####--------------------------------------------------------------------------------------------------------------------------------------------
## This script runs the simulations across all scenarios: implemented using the computational shared facility (CSF) at University of Manchester
####--------------------------------------------------------------------------------------------------------------------------------------------

#Load the simulation functions
source("./00_Simulation_Functions.R") #for CSF
# source(here::here("Scripts", "00_Simulation_Functions.R")) #for local run
library(tidyverse)

# Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  P = 10,
  RhoX = c(0, 0.5),
  Y_prev = c(0.2, 0.5),
  R2_based_on_maxR2 = c(TRUE, FALSE),
  betas = list(c(rep(log(1.10), 6), rep(log(1.5), 2), rep(log(2), 2)),
               c(rep(log(1.10), 3), log(1.5), log(2), rep(0, 5)))
)

args <- commandArgs(trailingOnly = T) #pull in all arguments from the qsub file
s <- as.numeric(args[1]) #Note that we need to specify what class the argument is

# number of repeats per scenario
n_rep <- 500

set.seed(698137 + s)

simultion_results <- simulation_nruns_fnc(n_iter = n_rep,
                                          P = sims_parameters$P[s],
                                          RhoX = sims_parameters$RhoX[s],
                                          Y_prev = sims_parameters$Y_prev[s],
                                          beta_true = sims_parameters$betas[s][[1]],
                                          R2_based_on_maxR2 = sims_parameters$R2_based_on_maxR2[s])

simultion_results <- simultion_results %>% 
  mutate("Simulation_Scenario" = s,
         "P" = sims_parameters$P[s],
         "RhoX" = sims_parameters$RhoX[s],
         "Y_prev" = sims_parameters$Y_prev[s],
         "R2_based_on_maxR2" = sims_parameters$R2_based_on_maxR2[s],
         "beta_true" = ifelse(all(sims_parameters$betas[s][[1]] == c(rep(log(1.10), 6), rep(log(1.5), 2), rep(log(2), 2))),
                              "AllX",
                              "HalfX"),
         .after = "Model")

write_rds(simultion_results, path = paste("./simulation_results_", s, ".RDS", sep = ""))

warnings() #return any warnings from the simulation (saved within the CSF) - note: none are returned for above
