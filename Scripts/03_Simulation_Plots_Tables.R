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
## This script produces the plots and tables for the paper (simulation results)
####-----------------------------------------------------------------------------------------

#First run the summary script to produce the relevant summarises from the simulation results
source(here::here("Scripts", "02_Simulation_Summary.R"))


####-----------------------------------------------------------------------------------------
## Table: create the simulation scenario table
####-----------------------------------------------------------------------------------------

simulation_scenario_table <- sims_all %>%
  distinct(SimulationScenario, .keep_all = TRUE) %>%
  select(SimulationScenario, RhoX, Y_prev, R2_based_on_maxR2, Beta) %>%
  mutate(R2_based_on_maxR2 = factor(ifelse(R2_based_on_maxR2 == TRUE,
                                           "Sample Size based on max R2",
                                           "Sample Size based on population R2"),
                                    levels = c("Sample Size based on population R2",
                                               "Sample Size based on max R2")),
         Beta = fct_recode(Beta,
                           "All 10 predictors" = "allX",
                           "5 predictors and 5 noise terms" = "halfX")) %>%
  rename("Simulation Scenario" = "SimulationScenario",
         "Rho" = "RhoX",
         "Prevalence of Y" = "Y_prev",
         "R2 for Sample Size Calculation" = "R2_based_on_maxR2") 

####-----------------------------------------------------------------------------------------
## Table: Summary of required minimum sample sizes across simulation scenarios
####-----------------------------------------------------------------------------------------

Table_samplesize_summary_sims <- sims_all %>%
  group_by(SimulationScenario) %>%
  summarise_at("SampleSize",
               .funs = list("median" = median,
                            "mean" = mean,
                            "min" = min,
                            "max" = max)) %>%
  ungroup() %>%
  left_join(sims_all %>%
              group_by(SimulationScenario) %>%
              summarise_at("SampleSizeCriteria",
                           .funs = list(function(x) { #summarise the sample criteria, based on whether unique across scenarios and iterations
                             ifelse(length(unique(x)) == 1,
                                    paste("Uniquely criteria ", unique(x)),
                                    paste("Varied with min=", min(x), ", and max=", max(x)))
                           })) %>%
              ungroup(),
            by = "SimulationScenario") %>%
  rename("Simulation Scenario" = "SimulationScenario",
         "Sample Size Criteria" = "SampleSizeCriteria")


####-----------------------------------------------------------------------------------------
## Tables: Summary of mean (overall) performance across simulation scenarios
####-----------------------------------------------------------------------------------------

Tables_mean_summary <- mean_summary_sims %>%
  group_by(PerformanceMetric) %>%
  nest() %>%
  mutate("SummaryTable" = map(data, function(x) {
    x %>%
      mutate("Lower" = mean - (1.96*sd),
             "Upper" = mean + (1.96*sd),
             #Present mean and 95% CI ready for table output:
             "Value" = paste(round(mean, 3), 
                             " (", 
                             round(Lower, 3), 
                             ", ", 
                             round(Upper, 3), 
                             ")", sep = "")) %>%
      select(SimulationScenario, Model, Value) %>%
      pivot_wider(names_from = "Model",
                  values_from = "Value") %>%
      rename("Simulation Scenario" = "SimulationScenario"
             
             , "MLE" = "mle"
             , "Uniform closed-form" = "uniform"
             , "Uniform bootstrap" = "bootstrap"
             , "Firths" = "Firths"
             , "LASSO" = "lasso"
             , "Ridge" = "ridge"
             )
  }))



####-----------------------------------------------------------------------------------------
## Figure: distribution of sample size R2 vs. model R2
####-----------------------------------------------------------------------------------------

R2_true_model_plot <- sims_all %>%
  select(SimulationScenario, 
         R2_SampCalc,
         R2_mle:R2_ridge) %>%
  pivot_longer(cols = c(R2_mle:R2_ridge), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([R2]+)_([A-Za-z]+)") %>%
  select(-PerformanceMetric) %>%
  mutate(SimulationScenario = factor(paste("Simulation Scenario \n", SimulationScenario, sep = " "),
                                     levels = paste("Simulation Scenario \n", 1:max(sims_all$SimulationScenario), 
                                                    sep = " "))) %>%
  filter(Model == "uniform") %>%
  ggplot(aes(x = value, y = R2_SampCalc)) +
  facet_wrap(~SimulationScenario,
             ncol = 4, nrow = 4, as.table = TRUE) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  #coord_equal() +
  xlab(paste("Cox-Snell R2 of the model, upon validation")) + 
  ylab("R2 used for Sample Size calculation") +
  theme_bw(base_size = 12)
  
  


####-----------------------------------------------------------------------------------------
## Figures showing distribution of performance metrics across simulation scenarios
####-----------------------------------------------------------------------------------------

Box_violin_CITL_plot <- sims_all %>%
  select(SimulationScenario,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         starts_with("CITL_")) %>% 
  pivot_longer(cols = c(starts_with("CITL_")), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
  mutate(Model = fct_relevel(fct_recode(Model
                                        , "MLE" = "mle"
                                        , "Uniform closed-form" = "uniform"
                                        , "Uniform bootstrap" = "bootstrap"
                                        , "Firths" = "Firths"
                                        , "LASSO" = "lasso"
                                        , "Ridge" = "ridge"
                                        ),
                             "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
                             "LASSO", "Ridge"),
         
         SimulationScenario = factor(paste("Simulation Scenario", SimulationScenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$SimulationScenario), 
                                                    sep = " "))) %>%
  filter(Y_prev == 0.2,
         RhoX == 0) %>% #focus on a subset of simulation scenarios (similar results for P(y=1)=50% and rho=0.5)
  ggplot(aes(x = Model, y = value, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Calibration-in-the-large")



Box_violin_CalSlope_plot <- sims_all %>%
  select(SimulationScenario,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         starts_with("CalSlope_")) %>% 
  pivot_longer(cols = c(starts_with("CalSlope_")), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
  mutate(Model = fct_relevel(fct_recode(Model,
                                        "MLE" = "mle",
                                        "Uniform closed-form" = "uniform",
                                        "Uniform bootstrap" = "bootstrap",
                                        "Firths" = "Firths",
                                        "LASSO" = "lasso",
                                        "Ridge" = "ridge"),
                             "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
                             "LASSO", "Ridge"),
         
         SimulationScenario = factor(paste("Simulation Scenario", SimulationScenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$SimulationScenario), 
                                                    sep = " "))) %>%
  filter(Y_prev == 0.2,
         RhoX == 0) %>% #focus on a subset of simulation scenarios (similar results for P(y=1)=50% and rho=0.5)
  ggplot(aes(x = Model, y = value, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Calibration Slope")



Box_violin_AUC_plot <- sims_all %>%
  select(SimulationScenario,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         starts_with("AUC_")) %>% 
  pivot_longer(cols = c(starts_with("AUC_")), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
  mutate(Model = fct_relevel(fct_recode(Model,
                                        "MLE" = "mle",
                                        "Uniform closed-form" = "uniform",
                                        "Uniform bootstrap" = "bootstrap",
                                        "Firths" = "Firths",
                                        "LASSO" = "lasso",
                                        "Ridge" = "ridge"),
                             "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
                             "LASSO", "Ridge"),
         
         SimulationScenario = factor(paste("Simulation Scenario", SimulationScenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$SimulationScenario), 
                                                    sep = " "))) %>%
  filter(Y_prev == 0.2,
         RhoX == 0) %>% #focus on a subset of simulation scenarios (similar results for P(y=1)=50% and rho=0.5)
  ggplot(aes(x = Model, y = value, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("AUC")



Box_violin_R2_plot <- sims_all %>%
  select(SimulationScenario,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         R2_mle:R2_ridge) %>% 
  pivot_longer(cols = c(R2_mle:R2_ridge), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([R2]+)_([A-Za-z]+)") %>%
  mutate(Model = fct_relevel(fct_recode(Model,
                                        "MLE" = "mle",
                                        "Uniform closed-form" = "uniform",
                                        "Uniform bootstrap" = "bootstrap",
                                        "Firths" = "Firths",
                                        "LASSO" = "lasso",
                                        "Ridge" = "ridge"),
                             "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
                             "LASSO", "Ridge"),
         
         SimulationScenario = factor(paste("Simulation Scenario", SimulationScenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$SimulationScenario), 
                                                    sep = " "))) %>%
  filter(Y_prev == 0.2,
         RhoX == 0) %>% #focus on a subset of simulation scenarios (similar results for P(y=1)=50% and rho=0.5)
  ggplot(aes(x = Model, y = value, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Cox-Snell R2")


####-----------------------------------------------------------------------------------------
## Figure: distribution of shrinkage/penalisation estimates
####-----------------------------------------------------------------------------------------

Box_violin_SF_plot <- sims_all %>%
  select(SimulationScenario, 
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         starts_with("SF_")) %>% 
  pivot_longer(cols = c(starts_with("SF_")), 
               names_to = "metric_model", 
               values_to = "value") %>%
  extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
  select(-PerformanceMetric) %>%
  mutate(Model = fct_relevel(fct_recode(Model,
                                        "Uniform closed-form" = "uniform",
                                        "Uniform bootstrap" = "bootstrap",
                                        "LASSO" = "lasso",
                                        "Ridge" = "ridge"),
                             "Uniform closed-form", "Uniform bootstrap", 
                             "LASSO", "Ridge"),
         SimulationScenario = factor(SimulationScenario,
                                     levels = paste(1:max(sims_all$SimulationScenario)))) %>%
  #filter(Y_prev == 0.2) %>% #focus on a subset of simulation scenarios (similar results for P(y=1)=50% and rho=0.5)
  ggplot(aes(x = SimulationScenario, y = value, group = )) +
  facet_grid(Model~., scales = "free_y") +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  xlab("Simulation Scenario") +
  ylab("Shrinkage factor/ Penalty term")
