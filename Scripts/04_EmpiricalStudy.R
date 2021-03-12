# #######################################################################################################################

# Author of code: Glen P. Martin.

# This is code for a empirical study presented in a manuscript entitled: 
# Developing Clinical Prediction Models: the importance of penalization methods and quantifying bootstrap 
# variability even when adhering to minimum sample size recommendations
# Authors:
#   Glen P. Martin
#   Richard Riley
#   Gary S. Collins
#   Matthew Sperrin


# #######################################################################################################################

library(tidyverse)
library(furrr)
library(pmsampsize)
library(glmnet)
library(logistf)
library(pROC)


mimic_data <- read_csv(file = here::here("Data", "mimic_penalised_cpms_cohort.csv"), 
                       col_names = TRUE,
                       # Adjust data structures, as needed:
                       col_types = cols(
                         gender = col_factor(levels = c("M", "F")),
                         admission_type = col_factor(levels = c("ELECTIVE", "URGENT", "EMERGENCY")),
                         ethnicity_grouped = col_factor(levels = c("white", "black", "asian",
                                                                   "hispanic", "other", "unknown"))
                       ))


####----------------------------------------------------------------------------------------
## Data clean
####----------------------------------------------------------------------------------------

#Clean age - define a categorical version of age (not normally recommended but needed here 
# since anyone with an age over 89 has their age adjusted by mimic team, so only know >89
mimic_data <- mimic_data %>%
  mutate(age_grouped = factor(ifelse(age < 30, "<30",
                                     ifelse(age < 40, "<40",
                                            ifelse(age < 50, "<50",
                                                   ifelse(age < 60, "<60",
                                                          ifelse(age < 70, "<70",
                                                                 ifelse(age < 80, "<80", 
                                                                        ">80")))))),
                              levels = c("<30", "<40", "<50", "<60", "<70", "<80", ">80"))
  )


#####---------------------------------------------------------------------------------------
### Simplify some categorical variables prior to modelling
#####---------------------------------------------------------------------------------------

mimic_data <- mimic_data %>%
  mutate(admission_type = fct_recode(admission_type,
                                     "Emergency" = "EMERGENCY", 
                                     "Non-Emergency" = "ELECTIVE", 
                                     "Non-Emergency" = "URGENT"),
         ethnicity_grouped = fct_recode(ethnicity_grouped,
                                        "white" = "white", 
                                        "black" = "black", 
                                        "other" = "asian",
                                        "other" = "hispanic", 
                                        "other" = "other", 
                                        "unknown" = "unknown"))



####----------------------------------------------------------------------------------------
## Subset the variables and rows that we wish to use for modelling
####----------------------------------------------------------------------------------------

Analysis_Cohort <- mimic_data %>%
  select(subject_id, hadm_id, icustay_id,
         age_grouped,
         gender,
         admission_type,
         ethnicity_grouped,
         ends_with("_mean"),
         hospital_expire_flag) %>%
  # Consider complete case analysis on the mimic iii cohort, for simplicity of
  #   modelling and bootstrap:
  filter(complete.cases(.))



####----------------------------------------------------------------------------------------
## Create Table: Required Sample Size calculation  
####----------------------------------------------------------------------------------------
Y_prev <- mean(Analysis_Cohort$hospital_expire_flag)
MaxR2 <- 1 - exp((2*(((Y_prev*100)*log(Y_prev)) + ((100 - (Y_prev*100))*log(1 - (Y_prev))))) / 100)
R2_for_SampCalc <- 0.15 * MaxR2

P <- ncol(model.matrix(hospital_expire_flag ~ .,
                       data = Analysis_Cohort %>%
                         select(hospital_expire_flag,
                                age_grouped,
                                gender,
                                admission_type,
                                ethnicity_grouped,
                                ends_with("_mean"))))


MIMIC_Sample_Size_Calc_Table <- as.data.frame(pmsampsize(type = "b",
                                                         rsquared = R2_for_SampCalc,
                                                         parameters = P,
                                                         shrinkage = 0.9,
                                                         prevalence = Y_prev)$results_table) %>%
  select(Samp_size) %>%
  rownames_to_column("Criteria") %>%
  pivot_wider(names_from = "Criteria", values_from = "Samp_size") %>%
  mutate("Max_Criteria" = which.max(c(`Criteria 1`, `Criteria 2`, `Criteria 3`)))






####----------------------------------------------------------------------------------------
## Create a function that fits all penalisation methods (for application to the
#    mimic iii example only!). We can then call this for each bootstrap internal validation
####----------------------------------------------------------------------------------------
modelling_fnc <- function(Data_for_Model_Fit,
                          Data_for_Predictions) {
  # Inputs: Data_for_Model_Fit = a tibble that one wishes to
  #                               use to fit/develop all of the models
  #         Data_for_Predictions = a tibble that one wishes to 
  #                               apply the developed models to and
  #                               obtain predictions on (can be
  #                               same as development data which gives 
  #                               apprarent performance)
  
  if (is.list(Data_for_Model_Fit)) Data_for_Model_Fit <- bind_rows(Data_for_Model_Fit)
  if (is.list(Data_for_Predictions)) Data_for_Predictions <- bind_rows(Data_for_Predictions)
  
  ## Create a design matrix of the development and validation data 
  ##  (useful for some of the modelling steps below):
  Dev_DM <- model.matrix(hospital_expire_flag ~ .,
                         data = Data_for_Model_Fit %>%
                           select(hospital_expire_flag,
                                  age_grouped,
                                  gender,
                                  admission_type,
                                  ethnicity_grouped,
                                  ends_with("_mean")))
  Val_DM <- model.matrix(hospital_expire_flag ~ .,
                         data = Data_for_Predictions %>%
                           select(hospital_expire_flag,
                                  age_grouped,
                                  gender,
                                  admission_type,
                                  ethnicity_grouped,
                                  ends_with("_mean")))
  
  P <- ncol(Dev_DM) - 1 #number of predictor terms (minus 1 for intercept)
  
  #####-------------------------------FIT THE MODELS--------------------------------------------
  ## Unpenalisaed MLE
  MLE_CPM <- glm(hospital_expire_flag ~ .,
                 data = Data_for_Model_Fit %>%
                   select(hospital_expire_flag,
                          age_grouped,
                          gender,
                          admission_type,
                          ethnicity_grouped,
                          ends_with("_mean")),
                 family = binomial(link = "logit"))
  MLE_pi_insampl <- predict(MLE_CPM, 
                            newdata = Data_for_Model_Fit, 
                            type = "response") #in-sample predictions
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate(MLE_CPM_Pi = MLE_pi_insampl)
  MLE_pi_outsampl <- predict(MLE_CPM, 
                             newdata = Data_for_Predictions, 
                             type = "response") #out-of-sample predictions
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate(MLE_CPM_Pi = MLE_pi_outsampl)
  
  
  ## Closed-form uniform shrinkage
  LR_MLE <- -2 * (as.numeric(logLik(with(Data_for_Model_Fit,
                                         glm(hospital_expire_flag ~ 1,
                                             family = binomial(link = "logit"))))) - 
                    as.numeric(logLik(MLE_CPM)))
  uniform_shrinkage <- 1 + (P / 
                              (nrow(Data_for_Model_Fit) * log(1 - (1 - exp(-LR_MLE/nrow(Data_for_Model_Fit))))))
  shrunkMLE_LP <- as.numeric((Dev_DM[,-1]) %*% 
                               (uniform_shrinkage*coef(MLE_CPM)[-1]))
  shrunk.MLE.beta <- c(coef(glm(Data_for_Model_Fit$hospital_expire_flag ~ offset(shrunkMLE_LP), 
                                family = binomial(link = "logit")))[1],
                       (uniform_shrinkage*coef(MLE_CPM)[-1])) #re-estimate the intercept of shrunk model
  #Apply shrunk model to calculate predictions:
  shrunk_CPM_validation_LP <- as.numeric(Dev_DM %*% shrunk.MLE.beta) #in-sample predictions
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate(UniformShrunk_CPM_Pi = (exp(shrunk_CPM_validation_LP)/(1 + exp(shrunk_CPM_validation_LP))))
  shrunk_CPM_validation_LP <- as.numeric(Val_DM %*% shrunk.MLE.beta) #out-of-sample predictions
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate(UniformShrunk_CPM_Pi = (exp(shrunk_CPM_validation_LP)/(1 + exp(shrunk_CPM_validation_LP))))
  
  
  ## Bootstrap uniform shrinkage
  Bootstrap_shrinkage <- rep(NA, 500)
  for (boot in 1:length(Bootstrap_shrinkage)) {
    BootData <- Data_for_Model_Fit %>%
      sample_n(n(), replace = TRUE)
    
    BootMod <- glm(hospital_expire_flag ~ .,
                   data = BootData %>%
                     select(hospital_expire_flag,
                            age_grouped,
                            gender,
                            admission_type,
                            ethnicity_grouped,
                            ends_with("_mean")),
                   family = binomial(link = "logit"))
    
    #Apply the bootstrap model to the original development data:
    Raw_LP <- Dev_DM %*% coef(BootMod)
    Bootstrap_shrinkage[boot] <- coef(glm(Data_for_Model_Fit$hospital_expire_flag ~ Raw_LP, 
                                          family = binomial(link = "logit")))[2] 
    
  }
  Bootstrap_shrinkage <- mean(Bootstrap_shrinkage)
  shrunkMLE_LP <- as.numeric((Dev_DM[,-1]) %*% (Bootstrap_shrinkage*coef(MLE_CPM)[-1]))
  bootshrunk.MLE.beta <- c(coef(glm(Data_for_Model_Fit$hospital_expire_flag ~ offset(shrunkMLE_LP), 
                                    family = binomial(link = "logit")))[1],
                           (Bootstrap_shrinkage*coef(MLE_CPM)[-1])) #re-estimate the intercept of shrunk model
  #Apply shrunk model to the validation set:
  shrunk_CPM_validation_LP <- as.numeric(Dev_DM %*% bootshrunk.MLE.beta) #in-sample predictions
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate(BootstrapShrunk_CPM_Pi = (exp(shrunk_CPM_validation_LP)/(1 + exp(shrunk_CPM_validation_LP))))
  shrunk_CPM_validation_LP <- as.numeric(Val_DM %*% bootshrunk.MLE.beta) #out-of-sample predictions
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate(BootstrapShrunk_CPM_Pi = (exp(shrunk_CPM_validation_LP)/(1 + exp(shrunk_CPM_validation_LP))))
  
  
  ## Firth's Bias-reduced logistic regression
  Firths_CPM <- logistf(hospital_expire_flag ~ .,
                        data = Data_for_Model_Fit %>%
                          select(hospital_expire_flag,
                                 age_grouped,
                                 gender,
                                 admission_type,
                                 ethnicity_grouped,
                                 ends_with("_mean")),
                        firth = TRUE)
  #re-estimate the intercept ready for prediction:
  Firths_LP <- as.numeric((Dev_DM[,-1]) %*% (as.numeric(coef(Firths_CPM))[-1]))
  Firths.beta <- c(coef(glm(Data_for_Model_Fit$hospital_expire_flag ~ offset(Firths_LP), 
                            family = binomial(link = "logit")))[1],
                   (as.numeric(coef(Firths_CPM))[-1])) 
  #Make predictions:
  Firths_CPM_validation_LP <- as.numeric(Dev_DM %*% Firths.beta) #in-sample predictions
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate(Firths_CPM_Pi = (exp(Firths_CPM_validation_LP)/(1 + exp(Firths_CPM_validation_LP))))
  Firths_CPM_validation_LP <- as.numeric(Val_DM %*% Firths.beta) #out-of-sample predictions
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate(Firths_CPM_Pi = (exp(Firths_CPM_validation_LP)/(1 + exp(Firths_CPM_validation_LP))))
  
  
  ## LASSO
  LASSO_CPM <- cv.glmnet(x = Dev_DM[,-1],
                         y = Data_for_Model_Fit$hospital_expire_flag,
                         family = "binomial",
                         alpha = 1,
                         nfolds = 10)
  #in-sample predictions:
  LASSO_Predictions <- predict(LASSO_CPM,
                               newx = Dev_DM[,-1],
                               s = "lambda.min",
                               type = "response")
  LASSO_Predictions1SE <- predict(LASSO_CPM,
                                  newx = Dev_DM[,-1],
                                  s = "lambda.1se",
                                  type = "response")
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate("LASSO_CPM_Pi" = c(LASSO_Predictions),
           "LASSO_CPM1SE_Pi" = c(LASSO_Predictions1SE))
  #out-of-sample predictions:
  LASSO_Predictions <- predict(LASSO_CPM,
                               newx = Val_DM[,-1],
                               s = "lambda.min",
                               type = "response")
  LASSO_Predictions1SE <- predict(LASSO_CPM,
                                  newx = Val_DM[,-1],
                                  s = "lambda.1se",
                                  type = "response")
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate("LASSO_CPM_Pi" = c(LASSO_Predictions),
           "LASSO_CPM1SE_Pi" = c(LASSO_Predictions1SE))
  
  
  ## Ridge
  Ridge_CPM <- cv.glmnet(x = Dev_DM[,-1],
                         y = Data_for_Model_Fit$hospital_expire_flag,
                         family = "binomial",
                         alpha = 0,
                         nfolds = 10)
  #in-sample predictions:
  Ridge_Predictions <- predict(Ridge_CPM,
                               newx = Dev_DM[,-1],
                               s = "lambda.min",
                               type = "response")
  Ridge_Predictions1SE <- predict(Ridge_CPM,
                                  newx = Dev_DM[,-1],
                                  s = "lambda.1se",
                                  type = "response")
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate("Ridge_CPM_Pi" = c(Ridge_Predictions),
           "Ridge_CPM1SE_Pi" = c(Ridge_Predictions1SE))
  #out-of-sample predictions:
  Ridge_Predictions <- predict(Ridge_CPM,
                               newx = Val_DM[,-1],
                               s = "lambda.min",
                               type = "response")
  Ridge_Predictions1SE <- predict(Ridge_CPM,
                                  newx = Val_DM[,-1],
                                  s = "lambda.1se",
                                  type = "response")
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate("Ridge_CPM_Pi" = c(Ridge_Predictions),
           "Ridge_CPM1SE_Pi" = c(Ridge_Predictions1SE))
  
  
  ## Repeat Cross-fold validation LASSO and Ridge
  Error_mat_LASSO <- NULL
  Error_mat_Ridge <- NULL
  for(r in 1:100) { #100 repeats of 10-fold cross validation
    RepeatCV_LASSO_CPM <- cv.glmnet(x = Dev_DM[,-1],
                                    y = Data_for_Model_Fit$hospital_expire_flag,
                                    family = "binomial",
                                    alpha = 1,
                                    nfolds = 10)
    Error_mat_LASSO <- Error_mat_LASSO %>% bind_rows(tibble("lambda" = RepeatCV_LASSO_CPM$lambda, 
                                                            "cvm" = RepeatCV_LASSO_CPM$cvm,
                                                            "cvm.se" = RepeatCV_LASSO_CPM$cvsd))
    
    RepeatCV_Ridge_CPM <- cv.glmnet(x = Dev_DM[,-1],
                                    y = Data_for_Model_Fit$hospital_expire_flag,
                                    family = "binomial",
                                    alpha = 0,
                                    nfolds = 10)
    Error_mat_Ridge <- Error_mat_Ridge %>% bind_rows(tibble("lambda" = RepeatCV_Ridge_CPM$lambda, 
                                                            "cvm" = RepeatCV_Ridge_CPM$cvm,
                                                            "cvm.se" = RepeatCV_Ridge_CPM$cvsd))
  }
  Summary_Error_mat_LASSO <- Error_mat_LASSO %>%
    group_by(lambda) %>%
    summarise(CVM = mean(cvm),
              se.CVM = mean(cvm.se),
              .groups = "drop") %>%
    arrange(desc(lambda)) %>%
    ungroup()
  Summary_Error_mat_Ridge <- Error_mat_Ridge %>%
    group_by(lambda) %>%
    summarise(CVM = mean(cvm),
              se.CVM = mean(cvm.se),
              .groups = "drop") %>%
    arrange(desc(lambda)) %>%
    ungroup()
  RepeatCV_LASSO_S <- Summary_Error_mat_LASSO$lambda[which.min(Summary_Error_mat_LASSO$CVM)]
  #in-sample predictions:
  RepeatCV_LASSO_Predictions <- predict(glmnet(x = Dev_DM[,-1],
                                               y = Data_for_Model_Fit$hospital_expire_flag,
                                               family = "binomial",
                                               alpha = 1),
                                        newx = Dev_DM[,-1],
                                        s = RepeatCV_LASSO_S,
                                        type = "response")
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate("RepeatCVLASSO_CPM_Pi" = c(RepeatCV_LASSO_Predictions))
  #out-of-sample predictions:
  RepeatCV_LASSO_Predictions <- predict(glmnet(x = Dev_DM[,-1],
                                               y = Data_for_Model_Fit$hospital_expire_flag,
                                               family = "binomial",
                                               alpha = 1),
                                        newx = Val_DM[,-1],
                                        s = RepeatCV_LASSO_S,
                                        type = "response")
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate("RepeatCVLASSO_CPM_Pi" = c(RepeatCV_LASSO_Predictions))
  
  RepeatCV_Ridge_S <- Summary_Error_mat_Ridge$lambda[which.min(Summary_Error_mat_Ridge$CVM)]
  #in-sample predictions:
  RepeatCV_Ridge_Predictions <- predict(glmnet(x = Dev_DM[,-1],
                                               y = Data_for_Model_Fit$hospital_expire_flag,
                                               family = "binomial",
                                               alpha = 0),
                                        newx = Dev_DM[,-1],
                                        s = RepeatCV_Ridge_S,
                                        type = "response")
  Data_for_Model_Fit <- Data_for_Model_Fit %>%
    mutate("RepeatCVRidge_CPM_Pi" = c(RepeatCV_Ridge_Predictions))
  #out-of-sample predictions:
  RepeatCV_Ridge_Predictions <- predict(glmnet(x = Dev_DM[,-1],
                                               y = Data_for_Model_Fit$hospital_expire_flag,
                                               family = "binomial",
                                               alpha = 0),
                                        newx = Val_DM[,-1],
                                        s = RepeatCV_Ridge_S,
                                        type = "response")
  Data_for_Predictions <- Data_for_Predictions %>%
    mutate("RepeatCVRidge_CPM_Pi" = c(RepeatCV_Ridge_Predictions))
  
  
  
  #####-------------------------------TEST THE MODELS--------------------------------------------
  ##Internal function to test predictive performance:
  Performance_fnc <- function(PredictedRisk, Response) {
    LP <- log(PredictedRisk / (1 - PredictedRisk))
    
    CITL_mod <- glm(Response ~ offset(LP), family = binomial(link = "logit"))
    CITL <- as.numeric(coef(CITL_mod)[1])
    CITL_se <- sqrt(vcov(CITL_mod)[1,1])
    
    CalSlope_mod <- glm(Response ~ LP, family = binomial(link = "logit"))
    CalSlope <- as.numeric(coef(CalSlope_mod)[2])
    CalSlope_se <- sqrt(vcov(CalSlope_mod)[2,2])
    
    AUC <- roc(response = Response,
               predictor = PredictedRisk,
               direction = "<",
               levels = c(0,1),
               ci = TRUE)
    
    BrierScore <- 1/length(Response) * (sum((PredictedRisk - Response)^2))
    
    return(list("CITL" = CITL,
                "CITL_se" = CITL_se,
                "CalSlope" = CalSlope,
                "CalSlope_se" = CalSlope_se,
                "AUC" = as.numeric(AUC$auc),
                "AUC_se" = sqrt(var(AUC)),
                "BrierScore" = BrierScore))
  }
  
  InSamplePredictivePerformance <- map(list("MLE" = Data_for_Model_Fit$MLE_CPM_Pi
                                             , "Uniform" = Data_for_Model_Fit$UniformShrunk_CPM_Pi
                                             , "Bootstrap" = Data_for_Model_Fit$BootstrapShrunk_CPM_Pi
                                             , "Firth" = Data_for_Model_Fit$Firths_CPM_Pi
                                             , "LASSO" = Data_for_Model_Fit$LASSO_CPM_Pi
                                             , "Ridge" = Data_for_Model_Fit$Ridge_CPM_Pi
                                             , "LASSO1SE" = Data_for_Model_Fit$LASSO_CPM1SE_Pi
                                             , "Ridge1SE" = Data_for_Model_Fit$Ridge_CPM1SE_Pi
                                             , "RepeatCVLASSO" = Data_for_Model_Fit$RepeatCVLASSO_CPM_Pi
                                             , "RepeatCVRidge" = Data_for_Model_Fit$RepeatCVRidge_CPM_Pi
                                             ),
                                        Performance_fnc,
                                        Response = Data_for_Model_Fit$hospital_expire_flag) %>%
    map_df(as_tibble, .id = "Model") %>%
    rename_at(.vars = vars(-Model),
              ~paste0("InSample_", ., sep = ""))
  OutSamplePredictivePerformance <- map(list("MLE" = Data_for_Predictions$MLE_CPM_Pi
                                             , "Uniform" = Data_for_Predictions$UniformShrunk_CPM_Pi
                                             , "Bootstrap" = Data_for_Predictions$BootstrapShrunk_CPM_Pi
                                             , "Firth" = Data_for_Predictions$Firths_CPM_Pi
                                             , "LASSO" = Data_for_Predictions$LASSO_CPM_Pi
                                             , "Ridge" = Data_for_Predictions$Ridge_CPM_Pi
                                             , "LASSO1SE" = Data_for_Predictions$LASSO_CPM1SE_Pi
                                             , "Ridge1SE" = Data_for_Predictions$Ridge_CPM1SE_Pi
                                             , "RepeatCVLASSO" = Data_for_Predictions$RepeatCVLASSO_CPM_Pi
                                             , "RepeatCVRidge" = Data_for_Predictions$RepeatCVRidge_CPM_Pi
                                             ),
                                        Performance_fnc,
                                        Response = Data_for_Predictions$hospital_expire_flag) %>%
    map_df(as_tibble, .id = "Model") %>%
    rename_at(.vars = vars(-Model),
              ~paste0("OutSample_", ., sep = ""))
  
  SF <- tibble("Model" = c("MLE", "Uniform", "Bootstrap", "Firth",
                           "LASSO", "Ridge", 
                           "LASSO1SE", "Ridge1SE", 
                           "RepeatCVLASSO", "RepeatCVRidge"),
               
               "SF" = c(NA, 
                        uniform_shrinkage, 
                        Bootstrap_shrinkage,
                        NA,
                        LASSO_CPM$lambda.min,
                        Ridge_CPM$lambda.min,
                        LASSO_CPM$lambda.1se,
                        Ridge_CPM$lambda.1se,
                        RepeatCV_LASSO_S,
                        RepeatCV_Ridge_S))
  
  #Save the performance results in a tibble:
  outputs <- InSamplePredictivePerformance %>%
    left_join(OutSamplePredictivePerformance, by = "Model") %>%
    left_join(SF, by = "Model")
  
  
  #####-------------------------------RETURN PREDICTED PERFORMANCE-------------------------------
  return(outputs)
}



####----------------------------------------------------------------------------------------
## Fit the models to the raw data, bootstrap data, and obtain internal validation results
####----------------------------------------------------------------------------------------


#Create a nested dataframe with each bootstrap data per (nested) row, 
# then pass this list of data to the above function
set.seed(260626)
nested_analysis_data <- tibble("Bootstrap_Index" = 0,
                               "Original_Data" = list(Analysis_Cohort),
                               #set first row to raw data (to obtain apparent performance):
                               "Bootstrap_Data" = list(Analysis_Cohort) 
                               ) %>%
  bind_rows(tibble("Bootstrap_Index" = 1:100, #set how many bootstrap sample we wish to take
                   "Original_Data" = list(Analysis_Cohort)) %>%
              #Apply the boostrap resampling:
              mutate("Bootstrap_Data" =  map(Original_Data,
                                             function(df) df %>% sample_n(nrow(df), replace = TRUE))))
 
#identical(nested_analysis_data$Bootstrap_Data[[20]], nested_analysis_data$Bootstrap_Data[[12]]) #=FALSE
#identical(nested_analysis_data$Bootstrap_Data[[10]], nested_analysis_data$Bootstrap_Data[[95]]) #=FALSE
#identical(nested_analysis_data$Original_Data[[10]], nested_analysis_data$Bootstrap_Data[[10]]) #=FALSE
#identical(nested_analysis_data$Original_Data[[1]], nested_analysis_data$Bootstrap_Data[[1]]) #=TRUE, as expected


##Apply the modelling functions to each bootstrap dataset (runs in parallel):
plan(multiprocess, workers = (availableCores() - 1))

nested_analysis_data <- nested_analysis_data %>%
  #Calculate the performance results of a model developed on each bootstrap data and tested in the original data:
  mutate("Results" = future_pmap(list(Data_for_Model_Fit = Bootstrap_Data,# fit models on bootstrap data
                                      Data_for_Predictions = Original_Data),#test them on the raw data (to get optimism)
                                 modelling_fnc,
                                 .progress = TRUE
                                 ) 
         )


# write_rds(nested_analysis_data, path = here::here("Data", "mimic_analysis_results.RDS"))



####----------------------------------------------------------------------------------------
## LOWER SAMPLE EXAMPLE: Fit the models to a sample of the raw data
####----------------------------------------------------------------------------------------


#Create a nested dataframe with each bootstrap data per (nested) row, 
# then pass this list of data to the above function
set.seed(98023)
Analysis_Cohort_Subset <- Analysis_Cohort %>%
  sample_n(MIMIC_Sample_Size_Calc_Table$Final, replace = FALSE)

subset_nested_analysis_data <- tibble("Bootstrap_Index" = 0,
                                      "Original_Data" = list(Analysis_Cohort_Subset),
                                      #set first row to raw data (to obtain apparent performance):
                                      "Bootstrap_Data" = list(Analysis_Cohort_Subset) 
                                      ) %>%
  bind_rows(tibble("Bootstrap_Index" = 1:100, #set how many bootstrap sample we wish to take
                   "Original_Data" = list(Analysis_Cohort_Subset)) %>%
              #Apply the boostrap resampling:
              mutate("Bootstrap_Data" =  map(Original_Data,
                                             function(df) df %>% sample_n(nrow(df), replace = TRUE))))

#identical(subset_nested_analysis_data$Bootstrap_Data[[20]], subset_nested_analysis_data$Bootstrap_Data[[12]]) #=FALSE
#identical(subset_nested_analysis_data$Bootstrap_Data[[33]], subset_nested_analysis_data$Bootstrap_Data[[75]]) #=FALSE
#identical(subset_nested_analysis_data$Original_Data[[44]], subset_nested_analysis_data$Bootstrap_Data[[44]]) #=FALSE
#identical(subset_nested_analysis_data$Original_Data[[1]], subset_nested_analysis_data$Bootstrap_Data[[1]]) #=TRUE, as expected

##Apply the modelling functions to each bootstrap dataset (runs in parallel):
plan(multiprocess, workers = (availableCores() - 1))

subset_nested_analysis_data <- subset_nested_analysis_data %>%
  #Calculate the performance results of a model developed on each bootstrap data and tested in the original data:
  mutate("Results" = future_pmap(list(Data_for_Model_Fit = Bootstrap_Data,# fit models on bootstrap data
                                      Data_for_Predictions = Original_Data),#test them on the raw data (to get optimism)
                                 modelling_fnc,
                                 .progress = TRUE)
         )


# write_rds(subset_nested_analysis_data, path = here::here("Data", "mimic_analysis_subset_results.RDS"))




####----------------------------------------------------------------------------------------
## Summarise results for entry into the manuscript
####----------------------------------------------------------------------------------------

# nested_analysis_data <- read_rds(here::here("Data", "mimic_analysis_results.RDS"))
# subset_nested_analysis_data <- read_rds(here::here("Data", "mimic_analysis_subset_results.RDS"))

# Unnest each result to create a tibble of performance results across the bootstraps:
summary_mimic_results <- nested_analysis_data %>%
  select(Bootstrap_Index, Results) %>%
  unnest(cols = c(Results)) %>%
  mutate("Study" = "Main", .before = "Bootstrap_Index") %>%
  bind_rows(subset_nested_analysis_data %>%
              select(Bootstrap_Index, Results) %>%
              unnest(cols = c(Results)) %>%
              mutate("Study" = "Subset", .before = "Bootstrap_Index"))

### Show the distribution in penalisation/ SF estimates across bootstrap samples
SF_Data <- summary_mimic_results %>%
  select(Study,
         Bootstrap_Index,
         Model,
         SF) %>% 
  filter(Bootstrap_Index != 0) %>%
  filter(Model != "MLE",
         Model != "Firth") %>%
  mutate(Model = fct_relevel(fct_recode(Model
                                        , "Uniform closed-form" = "Uniform"
                                        , "Uniform bootstrap" = "Bootstrap"
                                        , "LASSO" = "LASSO"
                                        , "LASSO (1SE)" = "LASSO1SE"
                                        , "Repeat \n CV LASSO" = "RepeatCVLASSO"
                                        , "Ridge" = "Ridge"
                                        , "Ridge (1SE)" = "Ridge1SE"
                                        , "Repeat \n CV Ridge" = "RepeatCVRidge"
                                        ),
                             "Uniform closed-form", "Uniform bootstrap",
                             "LASSO", "LASSO (1SE)", "Repeat \n CV LASSO", 
                             "Ridge", "Ridge (1SE)", "Repeat \n CV Ridge"))
MIMIC_sf_PLOT <- cowplot::plot_grid(
  SF_Data %>%
    filter(Model == "Uniform closed-form" |
             Model == "Uniform bootstrap") %>%
    ggplot(aes(x = Model, y = SF, color = Model)) +
    facet_grid(~Study, scales = "fixed", labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    xlab("Model") +
    ylab("Shrinkage")
  
  ,
  
  SF_Data %>%
    filter(Model == "LASSO" |
             Model == "LASSO (1SE)" |
             Model == "Repeat \n CV LASSO") %>%
    ggplot(aes(x = Model, y = SF, color = Model)) +
    facet_grid(~Study, scales = "fixed", labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    xlab("Model") +
    ylab("Penalty \n parameter")
  
  ,
  
  SF_Data %>%
    filter(Model == "Ridge" |
             Model == "Ridge (1SE)" |
             Model == "Repeat \n CV Ridge") %>%
    ggplot(aes(x = Model, y = SF, color = Model)) +
    facet_grid(~Study, scales = "fixed", labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    xlab("Model") +
    ylab("Penalty \n parameter")
  
  , align = c("hv")
  , nrow = 3
  , labels = c("A", "B", "C")
)




##First manipulate the data in a format we can work with:
Bootstrap_InternalValidation <- summary_mimic_results %>%
  filter(Bootstrap_Index != 0) %>%
  select(Study,
         Bootstrap_Index,
         Model,
         InSample_CITL,
         InSample_CalSlope,
         InSample_AUC,
         InSample_BrierScore,
         OutSample_CITL,
         OutSample_CalSlope,
         OutSample_AUC,
         OutSample_BrierScore) %>%
  #pivot the bootstrap internal validation results into long format:
  pivot_longer(col = c(-Study, -Bootstrap_Index, -Model),
               names_to = "Sample_Metric") %>%
  separate(Sample_Metric, into = c("Sample", "Metric"), sep = "_") %>%
  #Place the Outsample bootstrap performance estimates 'next to' the in-bootstrap-sample estimates:
  pivot_wider(id_cols = c("Study", "Bootstrap_Index", "Model", "Metric"),
              names_from = "Sample",
              values_from = "value") %>%
  #Calculate the optimism for each performance metric, as being the difference between
  #the performance of each bootstrap model within the bootstrap data (insample performance)
  #and the performance of each bootstrap model within the raw data (outsample performance):
  mutate(optimism = (InSample - OutSample)) %>%
  select(Study,
         Bootstrap_Index,
         Model,
         Metric,
         optimism) %>%
  #Join the bootstrap internal validation results above, with the apparent results:
  left_join(summary_mimic_results %>%
              filter(Bootstrap_Index == 0) %>%
              select(Study,
                     Model,
                     #need only select these results since OutSample_. are identical by definition
                     InSample_CITL,
                     InSample_CalSlope,
                     InSample_AUC,
                     InSample_BrierScore) %>%
              #pivot this into long format to match the bootstrap resutls above:
              pivot_longer(cols = c(-Study, -Model),
                           names_to = "Sample_Metric",
                           values_to = "ApparentPerformance") %>%
              separate(Sample_Metric, into = c("Sample", "Metric"), sep = "_") %>%
              select(-Sample),
            by = c("Study", "Model", "Metric")) %>%
  #Calculate the adjusted performance as the performance of the models fit
  #to the raw (original) data (i.e. apparent performance) minus each
  #optimism estimate across bootstraps:
  mutate(AdjustedPerformance = ApparentPerformance - optimism)

##Examine the mean performance (with 95% CI of this mean) across
#the bootstrap internal validation (i.e. average performance):
Tables_mimiciii_mean_summary <- Bootstrap_InternalValidation %>%
  mutate(Model = fct_relevel(fct_recode(Model
                                        , "MLE" = "MLE"
                                        , "Uniform closed-form" = "Uniform"
                                        , "Uniform bootstrap" = "Bootstrap"
                                        , "Firths" = "Firth"
                                        , "LASSO" = "LASSO"
                                        , "Repeat CV LASSO" = "RepeatCVLASSO"
                                        , "Ridge" = "Ridge"
                                        , "Repeat CV Ridge" = "RepeatCVRidge"
                                        ),
                             "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
                             "LASSO", "Repeat CV LASSO", "Ridge", "Repeat CV Ridge")) %>%
  group_by(Study, Model, Metric) %>%
  summarise("Mean" = mean(AdjustedPerformance),
            "sd" = sd(AdjustedPerformance),
            .groups = "drop") %>%
  ungroup() %>%
  mutate("Lower" = Mean - (1.96*sd),
         "Upper" = Mean + (1.96*sd),
         #Present mean and 95% CI ready for table output:
         "Value" = paste(round(Mean, 2),
                         " (",
                         round(Lower, 2),
                         ", ",
                         round(Upper, 2),
                         ")", sep = "")) %>%
  select(Study, Model, Metric, Value) %>%
  pivot_wider(id_cols = c("Study", "Model"),
              names_from = "Metric",
              values_from = "Value") %>%
  select(Study, Model, CITL, CalSlope, AUC, BrierScore) %>%
  rename("Calibration-in-the-large" = "CITL",
         "Calibration Slope" = "CalSlope",
         "AUC" = "AUC",
         "Brier Score" = "BrierScore")

#....we suggest that it is also constructive to look at the distribution of each
#adjusted performance results across the bootstrap samples, as well as the mean:
Box_violin_mimiciii_plot <- Bootstrap_InternalValidation %>%
  mutate(Model = fct_relevel(fct_recode(Model
                                        , "MLE" = "MLE"
                                        , "Uniform closed-form" = "Uniform"
                                        , "Uniform bootstrap" = "Bootstrap"
                                        , "Firths" = "Firth"
                                        , "LASSO" = "LASSO"
                                        , "Repeat CV LASSO" = "RepeatCVLASSO"
                                        , "Ridge" = "Ridge"
                                        , "Repeat CV Ridge" = "RepeatCVRidge"
                                        ),
                             "MLE", "Uniform closed-form", "Uniform bootstrap", "Firths", 
                             "LASSO", "Repeat CV LASSO", "Ridge", "Repeat CV Ridge"),

         Metric = fct_relevel(fct_recode(Metric,
                                         "Calibration-in-the-large" = "CITL",
                                         "Calibration Slope" = "CalSlope",
                                         "AUC" = "AUC",
                                         "Brier Score" = "BrierScore"),
                              "Calibration-in-the-large", "Calibration Slope", "AUC", "Brier Score")) %>%
  filter(Metric != "Brier Score") %>%
  ggplot(aes(x = Model, y = AdjustedPerformance, fill = Model)) +
  facet_wrap(~Metric + Study,
             nrow = 4, ncol = 2, scale = "free_y") +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1, color = "black") +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.1) +
  geom_blank(data = data.frame("Metric" = "Calibration-in-the-large",
                               "Model" = c("MLE", "Uniform closed-form", "Uniform bootstrap", "Firths",
                                           "LASSO", "Repeat CV LASSO", "Ridge", "Repeat CV Ridge"),
                               "AdjustedPerformance" = range(Bootstrap_InternalValidation$AdjustedPerformance[
                                 which(Bootstrap_InternalValidation$Metric=="CITL")
                               ]))) +
  geom_blank(data = data.frame("Metric" = "Calibration Slope",
                               "Model" = c("MLE", "Uniform closed-form", "Uniform bootstrap", "Firths",
                                           "LASSO", "Repeat CV LASSO", "Ridge", "Repeat CV Ridge"),
                               "AdjustedPerformance" = range(Bootstrap_InternalValidation$AdjustedPerformance[
                                 which(Bootstrap_InternalValidation$Metric=="CalSlope")
                               ]))) +
  geom_blank(data = data.frame("Metric" = "AUC",
                               "Model" = c("MLE", "Uniform closed-form", "Uniform bootstrap", "Firths",
                                           "LASSO", "Repeat CV LASSO", "Ridge", "Repeat CV Ridge"),
                               "AdjustedPerformance" = range(Bootstrap_InternalValidation$AdjustedPerformance[
                                 which(Bootstrap_InternalValidation$Metric=="AUC")
                               ]))) +
  # geom_blank(data = data.frame("Metric" = "Brier Score",
  #                              "Model" = c("MLE", "Uniform closed-form", "Uniform bootstrap", "Firths",
  #                                          "LASSO", "Repeat CV LASSO", "Ridge", "Repeat CV Ridge"),
  #                              "AdjustedPerformance" = range(Bootstrap_InternalValidation$AdjustedPerformance[
  #                                which(Bootstrap_InternalValidation$Metric=="BrierScore")
  #                              ]))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Bootstrap Adjusted Performance") +
  xlab("Estimation Method")


## Save the results for entry into the paper:
mimic_tables_figures <- list("MIMIC_Sample_Size_Calc_Table" = MIMIC_Sample_Size_Calc_Table,
                             "Tables_mimiciii_mean_summary" = Tables_mimiciii_mean_summary,
                             "MIMIC_sf_PLOT" = MIMIC_sf_PLOT,
                             "Box_violin_mimiciii_plot" = Box_violin_mimiciii_plot)

write_rds(mimic_tables_figures, path = here::here("Data", "mimic_tables_figures.RDS"))
