# #######################################################################################################################

# Author of code: Glen P. Martin.

# This is code for a simulation study presented in a manuscript entitled: 
# Developing Clinical Prediction Models: the importance of penalization methods and quantifying bootstrap 
# variability even when adhering to minimum sample size recommendations
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
  select(Simulation_Scenario, RhoX, Y_prev, R2_based_on_maxR2, Beta) %>%
  distinct(Simulation_Scenario, .keep_all = TRUE) %>%
  filter(RhoX == 0) %>%
  select(-RhoX) %>%
  mutate(R2_based_on_maxR2 = factor(ifelse(R2_based_on_maxR2 == TRUE,
                                           "Sample Size based on max R2",
                                           "Sample Size based on population R2"),
                                    levels = c("Sample Size based on population R2",
                                               "Sample Size based on max R2")),
         Beta = fct_recode(Beta,
                           "All 10 predictors" = "AllX",
                           "5 predictors and 5 noise terms" = "HalfX")) %>%
  rename("Simulation Scenario" = "Simulation_Scenario",
         "Prevalence of Y" = "Y_prev",
         "R2 for Sample Size Calculation" = "R2_based_on_maxR2") 

####-----------------------------------------------------------------------------------------
## Table: Summary of required minimum sample sizes across simulation scenarios
####-----------------------------------------------------------------------------------------

Table_samplesize_summary_sims <- sims_all %>%
  filter(RhoX == 0) %>%
  group_by(Simulation_Scenario) %>%
  summarise_at("SampleSize",
               .funs = list("median" = median,
                            "mean" = mean,
                            "min" = min,
                            "max" = max)) %>%
  ungroup() %>%
  left_join(sims_all %>%
              group_by(Simulation_Scenario) %>%
              summarise_at("SampleSizeCriteria",
                           .funs = list(function(x) { #summarise the sample criteria, based on whether unique across scenarios and iterations
                             ifelse(length(unique(x)) == 1,
                                    paste("Uniquely criteria ", unique(x)),
                                    paste("Varied with min=", min(x), ", and max=", max(x)))
                           })) %>%
              ungroup(),
            by = "Simulation_Scenario") %>%
  rename("Simulation Scenario" = "Simulation_Scenario",
         "Sample Size Criteria" = "SampleSizeCriteria")


####-----------------------------------------------------------------------------------------
## Tables: Summary of median (overall) performance across simulation scenarios
####-----------------------------------------------------------------------------------------

Tables_mean_summary <- mean_summary_sims %>%
  filter(Simulation_Scenario %in% 1:8) %>%
  group_by(PerformanceMetric) %>%
  nest() %>%
  mutate("SummaryTable" = map(data, function(x) {
    x %>%
      mutate(#Present mean and 95% CI ready for table output:
             "Value" = paste(format(round(median, 2), nsmall = 2),
                             " (", 
                             format(round(LQ, 2), nsmall = 2),
                             ", ", 
                             format(round(UQ, 2), nsmall = 2),
                             ")", sep = "")) %>%
      select(Simulation_Scenario, Model, Value) %>%
      pivot_wider(names_from = "Model",
                  values_from = "Value") %>%
      rename("Simulation Scenario" = "Simulation_Scenario")
    }
    )
    )



####-----------------------------------------------------------------------------------------
## Figure: distribution of sample size R2 vs. model R2
####-----------------------------------------------------------------------------------------

R2_true_model_plot <- sims_all %>%
  select(Simulation_Scenario, 
         Model,
         R2_SampCalc,
         R2) %>%
  filter(Simulation_Scenario %in% 1:8) %>%
  mutate(SimulationScenario = factor(paste("Simulation Scenario \n", Simulation_Scenario, sep = " "),
                                     levels = paste("Simulation Scenario \n", 1:max(sims_all$Simulation_Scenario), 
                                                    sep = " "))) %>%
  filter(Model == "MLE") %>%
  ggplot(aes(x = R2, y = R2_SampCalc)) +
  facet_wrap(~SimulationScenario,
             ncol = 4, nrow = 4, as.table = TRUE) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  #coord_equal() +
  xlab(paste("Cox-Snell R2 of the model, upon validation")) + 
  ylab("R2 used for Sample Size calculation") +
  theme_bw(base_size = 12)
  
  


####-----------------------------------------------------------------------------------------
## Figures showing distribution of performance metrics across simulation scenarios
####-----------------------------------------------------------------------------------------

Box_violin_CITL_plot <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         CITL) %>% 
  mutate(Diff = (CITL - 0)^2,
         SimulationScenario = factor(paste("Simulation Scenario", Simulation_Scenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$Simulation_Scenario), 
                                                    sep = " "))) %>%
  group_by(Model, Simulation_Scenario) %>%
  mutate(RMSD = format( round(sqrt(mean(Diff)),2) , nsmall = 2)) %>%
  ungroup() %>%
  filter(RhoX == 0) %>% 
  ggplot(aes(x = Model, y = CITL, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_text(aes(y = max(sims_all$CITL), label=RMSD), color = "black", fontface = 2, size = 3) +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Calibration-in-the-large")



Box_violin_CalSlope_plot <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         CalSlope) %>% 
  mutate(Diff = (CalSlope - 1)^2,
         SimulationScenario = factor(paste("Simulation Scenario", Simulation_Scenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$Simulation_Scenario), 
                                                    sep = " "))) %>%
  group_by(Model, Simulation_Scenario) %>%
  mutate(RMSD = format( round(sqrt(mean(Diff)),2) , nsmall = 2)) %>%
  ungroup() %>%
  filter(RhoX == 0) %>% 
  ggplot(aes(x = Model, y = CalSlope, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_text(aes(y = max(sims_all$CalSlope), label=RMSD), color = "black", fontface = 2, size = 3) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) + 
  geom_jitter(shape = 16, position = position_jitter(0.4), alpha = 0.1, aes(colour = Model)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0.6, 1.8, by = 0.2)) +
  ylab("Calibration Slope")



Box_violin_AUC_plot <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         AUC) %>% 
  mutate(SimulationScenario = factor(paste("Simulation Scenario", Simulation_Scenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$Simulation_Scenario), 
                                                    sep = " "))) %>%
  filter(RhoX == 0) %>% 
  ggplot(aes(x = Model, y = AUC, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("AUC")



Box_violin_R2_plot <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         R2) %>% 
  mutate(SimulationScenario = factor(paste("Simulation Scenario", Simulation_Scenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$Simulation_Scenario), 
                                                    sep = " "))) %>%
  filter(RhoX == 0) %>%  
  ggplot(aes(x = Model, y = R2, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Cox-Snell R2")



Box_violin_BrierScore_plot <- sims_all %>%
  select(Simulation_Scenario,
         Model,
         P, RhoX, Y_prev, R2_based_on_maxR2, Beta,
         BrierScore) %>% 
  mutate(SimulationScenario = factor(paste("Simulation Scenario", Simulation_Scenario, sep = " "),
                                     levels = paste("Simulation Scenario", 1:max(sims_all$Simulation_Scenario), 
                                                    sep = " "))) %>%
  filter(RhoX == 0) %>% 
  ggplot(aes(x = Model, y = BrierScore, fill = Model)) +
  facet_wrap(~SimulationScenario, scales = "fixed", 
             ncol = 2, as.table = TRUE) +
  geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Brier Score")







####-----------------------------------------------------------------------------------------
## Figure: distribution of shrinkage/penalisation estimates
####-----------------------------------------------------------------------------------------

Box_violin_SF_plot <- cowplot::plot_grid(
  sims_all %>%
    select(Simulation_Scenario, 
           Beta,
           SF_uniform, SF_bootstrap) %>% 
    pivot_longer(cols = c(starts_with("SF_")), 
                 names_to = "metric_model", 
                 values_to = "value") %>%
    extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z]+)") %>%
    select(-PerformanceMetric) %>%
    mutate(Model = fct_relevel(fct_recode(Model,
                                          "Uniform closed-form" = "uniform",
                                          "Uniform bootstrap" = "bootstrap"
                                          ),
                               "Uniform closed-form", "Uniform bootstrap"),
    SimulationScenario = factor(Simulation_Scenario,
                                levels = paste(1:max(sims_all$Simulation_Scenario)))) %>%
    filter(SimulationScenario %in% 1:8) %>%
    ggplot(aes(x = SimulationScenario, y = value, color = Beta)) +
    facet_grid(Model~., scales = "fixed", labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    xlab("Simulation Scenario") +
    ylab("Shrinkage")
  ,
  sims_all %>%
    select(Simulation_Scenario, 
           Beta,
           SF_lasso, SF_ReapeatCVlasso) %>% 
    pivot_longer(cols = c(starts_with("SF_")), 
                 names_to = "metric_model", 
                 values_to = "value") %>%
    extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z1]+)") %>%
    select(-PerformanceMetric) %>%
    mutate(Model = fct_relevel(fct_recode(Model
                                          , "LASSO" = "lasso"
                                          , "Repeat CV LASSO" = "ReapeatCVlasso"
                                          ),
                               "LASSO", "Repeat CV LASSO"),
    SimulationScenario = factor(Simulation_Scenario,
                                levels = paste(1:max(sims_all$Simulation_Scenario)))) %>%
    filter(SimulationScenario %in% 1:8) %>%
    ggplot(aes(x = SimulationScenario, y = value, color = Beta)) +
    facet_grid(Model~., scales = "fixed", labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    xlab("Simulation Scenario") +
    ylab("Penalty parameter")
  ,
  sims_all %>%
    select(Simulation_Scenario, 
           Beta,
           SF_ridge, SF_ReapeatCVridge) %>% 
    pivot_longer(cols = c(starts_with("SF_")), 
                 names_to = "metric_model", 
                 values_to = "value") %>%
    extract(metric_model, into = c("PerformanceMetric", "Model"), "([A-Za-z]+)_([A-Za-z1]+)") %>%
    select(-PerformanceMetric) %>%
    mutate(Model = fct_relevel(fct_recode(Model
                                          , "Ridge" = "ridge"
                                          , "Repeat CV Ridge" = "ReapeatCVridge"
                                          ),
                               "Ridge", "Repeat CV Ridge"),
    SimulationScenario = factor(Simulation_Scenario,
                                levels = paste(1:max(sims_all$Simulation_Scenario)))) %>%
    filter(SimulationScenario %in% 1:8) %>%
    ggplot(aes(x = SimulationScenario, y = value, color = Beta)) +
    facet_grid(Model~., scales = "fixed", labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    xlab("Simulation Scenario") +
    ylab("Penalty parameter")
  
  , align = c("hv")
  , nrow = 3
  , labels = c("A", "B", "C")
)

