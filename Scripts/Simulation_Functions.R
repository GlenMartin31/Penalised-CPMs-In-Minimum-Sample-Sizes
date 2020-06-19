####-----------------------------------------------------------------------------------------
## Define a function to run the proccesses within a single iteration
####-----------------------------------------------------------------------------------------
simulation_singlerun_fnc <- function(P,
                                     Sigma.mat,
                                     Y_prev,
                                     beta_true,
                                     R2_based_on_maxR2) {
  #Input: P = the number of predictors to simulate 
  #       Sigma.mat = covariance matrix to generate X
  #       Y_prev = the overall event rate in the population
  #       beta_true = a vector of 'true' log odds ratios to generate the binary outcomes
  #       R2_based_on_maxR2 = indicator of whether the sample size calculation is based on
  #             15% of the maximum R2 for given y_prev (TRUE), or based on population dervied R2 (FALSE)
  
  require(tidyverse)
  require(pmsampsize)
  require(glmnet)
  require(logistf)
  require(pROC)
  
  N_pop <- 1000000
  X <- MASS::mvrnorm(n = N_pop, mu = rep(0, P), Sigma = Sigma.mat)
  
  beta <- c(as.numeric(coef(glm(rbinom(N_pop, 1, prob = Y_prev) ~ 
                                  offset(as.vector(X %*% beta_true)), 
                                family = binomial(link = "logit")))[1]),
            beta_true) #set intercept to control overall event rate
  
  colnames(X) <- paste("V", 1:P, sep = "")
  PopulationData <- tibble(
    ID = 1:N_pop,
    as_tibble(X),
    LP = as.vector(cbind(1, X) %*% beta),
    Pi = (exp(LP)/(1 + exp(LP))),
    Y = rbinom(n = N_pop, 1, Pi)
    )
  
  ### Calculate R2 for the sample size calculation:
  if (R2_based_on_maxR2 == TRUE) {
    MaxR2 <- 1 - exp((2*(((Y_prev*100)*log(Y_prev)) + ((100 - (Y_prev*100))*log(1 - (Y_prev))))) / 100)
    R2_for_SampCalc <- 0.15 * MaxR2
  }else{
    #fit a model to the population data to estimate a value of R2 for the sample size calculation
    PopModel <- with(PopulationData,
                     glm(as.formula(paste("Y ~ ", paste((paste("V", 1:P, sep = "", collapse = " + "))))),
                         family = binomial(link = "logit")))
    LR <- -2 * (as.numeric(logLik(glm(Y ~ 1, data = PopulationData,
                                      family = binomial(link = "logit")))) - 
                  as.numeric(logLik(PopModel)))
    R2_app <- 1 - exp(-LR/N_pop)
    S_vh <- 1 + (length(coef(PopModel)) / (N_pop*log(1 - R2_app)))
    R2_for_SampCalc <- S_vh*R2_app
  }
  ###Calculate the required sample size
  SampleSize <- as.data.frame(pmsampsize(type = "b",
                                         rsquared = R2_for_SampCalc,
                                         parameters = P,
                                         shrinkage = 0.9,
                                         prevalence = Y_prev)$results_table) %>% 
                select(Samp_size) %>%
                rownames_to_column("Criteria") %>%
                pivot_wider(names_from = "Criteria", values_from = "Samp_size") %>%
                mutate("Max_Criteria" = which.max(c(`Criteria 1`, `Criteria 2`, `Criteria 3`)))
  N_Dev <- SampleSize$Final
  
  #Randomly select a derivation dataset from the population data
  DevelopmentData <- PopulationData %>%
    sample_n(size = N_Dev,
             replace = FALSE)
  ValidationData <- PopulationData %>%
    filter(ID %in% DevelopmentData$ID == FALSE) #use the remaining (unsampled) data as a validation set
  
  ### Fit each model to the development data
  ## Unpenalisaed MLE
  MLE_CPM <- with(DevelopmentData,
                  glm(as.formula(paste("Y ~ ", paste((paste("V", 1:P, sep = "", collapse = " + "))))),
                      family = binomial(link = "logit")))
  MLE_pi <- predict(MLE_CPM, newdata = ValidationData, type = "response")
  ValidationData <- ValidationData %>%
    mutate(MLE_CPM_Pi = MLE_pi)
  
  
  ## Closed-form uniform shrinkage
  LR_MLE <- -2 * (as.numeric(logLik(with(DevelopmentData,
                                         glm(Y ~ 1,
                                             family = binomial(link = "logit"))))) - 
                    as.numeric(logLik(MLE_CPM)))
  uniform_shrinkage <- 1 + (P / 
                              (nrow(DevelopmentData) * log(1 - (1 - exp(-LR_MLE/nrow(DevelopmentData))))))
  #re-estimate the intercept of shrunk model:
  shrunkMLE_LP <- as.numeric((DevelopmentData %>%
                                select(starts_with("V")) %>%
                                data.matrix()) %*% (uniform_shrinkage*coef(MLE_CPM)[-1]))
  shrunk.MLE.beta <- c(coef(glm(DevelopmentData$Y ~ offset(shrunkMLE_LP), 
                                family = binomial(link = "logit")))[1],
                       (uniform_shrinkage*coef(MLE_CPM)[-1])) 
  #Apply shrunk model to the validation set:
  shrunk_CPM_validation_LP <- as.numeric(cbind(1,
                                               ValidationData %>%
                                                 select(starts_with("V")) %>%
                                                 data.matrix()) %*% shrunk.MLE.beta) 
  ValidationData <- ValidationData %>%
    mutate(UniformShrunk_CPM_Pi = (exp(shrunk_CPM_validation_LP)/(1 + exp(shrunk_CPM_validation_LP))))
  
  
  ## Bootstrap uniform shrinkage
  Bootstrap_shrinkage <- rep(NA, 500)
  for (boot in 1:500) {
    BootData <- DevelopmentData %>%
      sample_n(n(), replace = TRUE)
    
    BootMod <- with(BootData,
                    glm(as.formula(paste("Y ~ ", paste((paste("V", 1:P, sep = "", collapse = " + "))))),
                        family = binomial(link = "logit")))
    
    Raw_LP <- cbind(1,
                    DevelopmentData %>%
                      select(starts_with("V")) %>%
                      data.matrix()) %*% coef(BootMod)
    Bootstrap_shrinkage[boot] <- 1 - #"1 minus" since BootMod perfectly calibrated, by definition, in BootData
      coef(glm(DevelopmentData$Y ~ Raw_LP, 
               family = binomial(link = "logit")))[2] 
    
  }
  Bootstrap_shrinkage <- 1 - mean(Bootstrap_shrinkage)
  #re-estimate the intercept of shrunk model:
  shrunkMLE_LP <- as.numeric((DevelopmentData %>%
                                select(starts_with("V")) %>%
                                data.matrix()) %*% (Bootstrap_shrinkage*coef(MLE_CPM)[-1]))
  bootshrunk.MLE.beta <- c(coef(glm(DevelopmentData$Y ~ offset(shrunkMLE_LP), 
                                    family = binomial(link = "logit")))[1],
                           (Bootstrap_shrinkage*coef(MLE_CPM)[-1])) 
  #Apply shrunk model to the validation set:
  shrunk_CPM_validation_LP <- as.numeric(cbind(1,
                                               ValidationData %>%
                                                 select(starts_with("V")) %>%
                                                 data.matrix()) %*% bootshrunk.MLE.beta) 
  ValidationData <- ValidationData %>%
    mutate(BootstrapShrunk_CPM_Pi = (exp(shrunk_CPM_validation_LP)/(1 + exp(shrunk_CPM_validation_LP))))
  
  
  ## Firth's Bias-reduced logistic regression
  Firths_CPM <- logistf(as.formula(paste("Y ~ ", paste((paste("V", 1:P, sep = "", collapse = " + "))))),
                        data = DevelopmentData,
                        firth = TRUE)
  Firths_CPM_validation_LP <- as.numeric(cbind(1,
                                               ValidationData %>%
                                                 select(starts_with("V")) %>%
                                                 data.matrix()) %*% as.numeric(coef(Firths_CPM)))
  ValidationData <- ValidationData %>%
    mutate(Firths_CPM_Pi = (exp(Firths_CPM_validation_LP)/(1 + exp(Firths_CPM_validation_LP))))
  
  
  ## LASSO
  LASSO_CPM <- cv.glmnet(x = DevelopmentData %>%
                           select(starts_with("V")) %>%
                           data.matrix(),
                         y = DevelopmentData$Y,
                         family = "binomial",
                         alpha = 1,
                         nfolds = 10)
  LASSO_Predictions <- predict(LASSO_CPM,
                               newx = ValidationData %>%
                                 select(starts_with("V")) %>%
                                 data.matrix(),
                               s = "lambda.min",
                               type = "response")
  ValidationData <- ValidationData %>%
    mutate("LASSO_CPM_Pi" = c(LASSO_Predictions))
  
  
  ## Ridge
  Ridge_CPM <- cv.glmnet(x = DevelopmentData %>%
                           select(starts_with("V")) %>%
                           data.matrix(),
                         y = DevelopmentData$Y,
                         family = "binomial",
                         alpha = 0,
                         nfolds = 10)
  Ridge_Predictions <- predict(Ridge_CPM,
                               newx = ValidationData %>%
                                 select(starts_with("V")) %>%
                                 data.matrix(),
                               s = "lambda.min",
                               type = "response")
  ValidationData <- ValidationData %>%
    mutate("Ridge_CPM_Pi" = c(Ridge_Predictions))
  
  
  ### Test each model in the validation data
  PredictivePerformance <- map(list("MLE" = ValidationData$MLE_CPM_Pi
                                    , "Uniform" = ValidationData$UniformShrunk_CPM_Pi
                                    , "Bootstrap" = ValidationData$BootstrapShrunk_CPM_Pi
                                    , "Firth" = ValidationData$Firths_CPM_Pi
                                    , "LASSO" = ValidationData$LASSO_CPM_Pi
                                    , "Ridge" = ValidationData$Ridge_CPM_Pi
                                    ),
                               Performance_fnc,
                               Response = ValidationData$Y)
  
  ### Extract the results
  outputs <- list()
  outputs$SampleSize <- N_Dev
  outputs$SampleSizeCriteria <- SampleSize$Max_Criteria
  
  outputs$SF_uniform <- uniform_shrinkage
  outputs$SF_bootstrap <- Bootstrap_shrinkage
  outputs$SF_lasso <- LASSO_CPM$lambda.min
  outputs$SF_ridge <- Ridge_CPM$lambda.min
  
  outputs$R2_SampCalc <- R2_for_SampCalc
  outputs$R2_mle <- PredictivePerformance$MLE$R2_coxsnell
  outputs$R2_uniform <- PredictivePerformance$Uniform$R2_coxsnell
  outputs$R2_bootstrap <- PredictivePerformance$Bootstrap$R2_coxsnell
  outputs$R2_Firths <- PredictivePerformance$Firth$R2_coxsnell
  outputs$R2_lasso <- PredictivePerformance$LASSO$R2_coxsnell
  outputs$R2_ridge <- PredictivePerformance$Ridge$R2_coxsnell
  
  outputs$CITL_mle <- PredictivePerformance$MLE$CITL
  outputs$CITLLower_mle <- PredictivePerformance$MLE$CITL_lower
  outputs$CITLUpper_mle <- PredictivePerformance$MLE$CITL_upper
  outputs$CITL_uniform <- PredictivePerformance$Uniform$CITL
  outputs$CITLLower_uniform <- PredictivePerformance$Uniform$CITL_lower
  outputs$CITLUpper_uniform <- PredictivePerformance$Uniform$CITL_upper
  outputs$CITL_bootstrap <- PredictivePerformance$Bootstrap$CITL
  outputs$CITLLower_bootstrap <- PredictivePerformance$Bootstrap$CITL_lower
  outputs$CITLUpper_bootstrap <- PredictivePerformance$Bootstrap$CITL_upper
  outputs$CITL_Firths <- PredictivePerformance$Firth$CITL
  outputs$CITLLower_Firths <- PredictivePerformance$Firth$CITL_lower
  outputs$CITLUpper_Firths <- PredictivePerformance$Firth$CITL_upper
  outputs$CITL_lasso <- PredictivePerformance$LASSO$CITL
  outputs$CITLLower_lasso <- PredictivePerformance$LASSO$CITL_lower
  outputs$CITLUpper_lasso <- PredictivePerformance$LASSO$CITL_upper
  outputs$CITL_ridge <- PredictivePerformance$Ridge$CITL
  outputs$CITLLower_ridge <- PredictivePerformance$Ridge$CITL_lower
  outputs$CITLUpper_ridge <- PredictivePerformance$Ridge$CITL_upper
  
  outputs$CalSlope_mle <- PredictivePerformance$MLE$CalSlope
  outputs$CalSlopeLower_mle <- PredictivePerformance$MLE$CalSlope_lower
  outputs$CalSlopeUpper_mle <- PredictivePerformance$MLE$CalSlope_upper
  outputs$CalSlope_uniform <- PredictivePerformance$Uniform$CalSlope
  outputs$CalSlopeLower_uniform <- PredictivePerformance$Uniform$CalSlope_lower
  outputs$CalSlopeUpper_uniform <- PredictivePerformance$Uniform$CalSlope_upper
  outputs$CalSlope_bootstrap <- PredictivePerformance$Bootstrap$CalSlope
  outputs$CalSlopeLower_bootstrap <- PredictivePerformance$Bootstrap$CalSlope_lower
  outputs$CalSlopeUpper_bootstrap <- PredictivePerformance$Bootstrap$CalSlope_upper
  outputs$CalSlope_Firths <- PredictivePerformance$Firth$CalSlope
  outputs$CalSlopeLower_Firths <- PredictivePerformance$Firth$CalSlope_lower
  outputs$CalSlopeUpper_Firths <- PredictivePerformance$Firth$CalSlope_upper
  outputs$CalSlope_lasso <- PredictivePerformance$LASSO$CalSlope
  outputs$CalSlopeLower_lasso <- PredictivePerformance$LASSO$CalSlope_lower
  outputs$CalSlopeUpper_lasso <- PredictivePerformance$LASSO$CalSlope_upper
  outputs$CalSlope_ridge <- PredictivePerformance$Ridge$CalSlope
  outputs$CalSlopeLower_ridge <- PredictivePerformance$Ridge$CalSlope_lower
  outputs$CalSlopeUpper_ridge <- PredictivePerformance$Ridge$CalSlope_upper
  
  outputs$AUC_mle <- PredictivePerformance$MLE$AUC
  outputs$AUCLower_mle <- PredictivePerformance$MLE$AUC_lower
  outputs$AUCUpper_mle <- PredictivePerformance$MLE$AUC_upper
  outputs$AUC_uniform <- PredictivePerformance$Uniform$AUC
  outputs$AUCLower_uniform <- PredictivePerformance$Uniform$AUC_lower
  outputs$AUCUpper_uniform <- PredictivePerformance$Uniform$AUC_upper
  outputs$AUC_bootstrap <- PredictivePerformance$Bootstrap$AUC
  outputs$AUCLower_bootstrap <- PredictivePerformance$Bootstrap$AUC_lower
  outputs$AUCUpper_bootstrap <- PredictivePerformance$Bootstrap$AUC_upper
  outputs$AUC_Firths <- PredictivePerformance$Firth$AUC
  outputs$AUCLower_Firths <- PredictivePerformance$Firth$AUC_lower
  outputs$AUCUpper_Firths <- PredictivePerformance$Firth$AUC_upper
  outputs$AUC_lasso <- PredictivePerformance$LASSO$AUC
  outputs$AUCLower_lasso <- PredictivePerformance$LASSO$AUC_lower
  outputs$AUCUpper_lasso <- PredictivePerformance$LASSO$AUC_upper
  outputs$AUC_ridge <- PredictivePerformance$Ridge$AUC
  outputs$AUCLower_ridge <- PredictivePerformance$Ridge$AUC_lower
  outputs$AUCUpper_ridge <- PredictivePerformance$Ridge$AUC_upper
  
  ### Return the results
  return(outputs)
}


####-----------------------------------------------------------------------------------------
## Define function to assess predictive performance of each model
####-----------------------------------------------------------------------------------------
Performance_fnc <- function(PredictedRisk, Response) {
  LP <- log(PredictedRisk / (1 - PredictedRisk))
  
  CITL_mod <- glm(Response ~ offset(LP), family = binomial(link = "logit"))
  CILT_confint <- confint(profile(CITL_mod), level = 0.95)
  CITL <- as.numeric(coef(CITL_mod)[1])
  CITL_lower <- CILT_confint["2.5 %"]
  CITL_upper <- CILT_confint["97.5 %"]
  
  CalSlope_mod <- glm(Response ~ LP, family = binomial(link = "logit"))
  CalSlope_confint <- confint(profile(CalSlope_mod), level = 0.95)
  CalSlope <- as.numeric(coef(CalSlope_mod)[2])
  CalSlope_lower <- as.numeric(CalSlope_confint["LP", "2.5 %"])
  CalSlope_upper <- as.numeric(CalSlope_confint["LP", "97.5 %"])
  
  AUC <- roc(response = Response,
             predictor = PredictedRisk,
             direction = "<",
             levels = c(0,1),
             ci = TRUE)
  
  LR <- -2 * (as.numeric(logLik(glm(Response ~ 1, family = binomial(link = "logit")))) - 
                as.numeric(logLik(CITL_mod)))
  R2_coxsnell <- 1 - exp(-LR/length(Response))
  
  return(list("CITL" = CITL,
              "CITL_lower" = CITL_lower,
              "CITL_upper" = CITL_upper,
              "CalSlope" = CalSlope,
              "CalSlope_lower" = CalSlope_lower,
              "CalSlope_upper" = CalSlope_upper,
              "AUC" = as.numeric(AUC$auc),
              "AUC_lower" = as.numeric(AUC$ci)[1],
              "AUC_upper" = as.numeric(AUC$ci)[3],
              "R2_coxsnell" = R2_coxsnell))
}



####-----------------------------------------------------------------------------------------
## Define a function to repeat the simulation across all iterations
####-----------------------------------------------------------------------------------------
simulation_nruns_fnc <- function(n_iter = 100,
                                 P,
                                 RhoX,
                                 Y_prev,
                                 beta_true,
                                 R2_based_on_maxR2) {
  #Input: n_iter = the number of iterations to repeat the simulation over
  #       RhoX = a vector of correlation values to run the simulation over: this is pairwise in that
  #                     each element of RhoX give cor(X1, X2), cor(X3, X4), etc.
  #       Y_prev = the overall event rate in the population
  #       beta_true = a vector of 'true' log odds ratios to generate the binary outcomes
  #       R2_based_on_maxR2 = indicator of whether the sample size calculation is based on
  #             15% of the maximum R2 for given y_prev (TRUE), or based on population dervied R2 (FALSE)
  
  require(tidyverse)
  
  Sigma.mat <- diag(1, ncol = P, nrow = P)
  #Create pair-wise correlations:
  if (P <= 2) {
    Sigma.mat[upper.tri(Sigma.mat)] <- Sigma.mat[lower.tri(Sigma.mat)] <- RhoX
  }else{
    Sigma.mat[(row(Sigma.mat) - col(Sigma.mat)) == 1 | 
                (row(Sigma.mat) - col(Sigma.mat)) == -1][c(TRUE, TRUE, FALSE, FALSE)] <- RhoX
  }
  
  ## Define a matrix to store the results across all iterations
  result_names <- c("SampleSize", "SampleSizeCriteria"
                    , "SF_uniform", "SF_bootstrap", "SF_lasso", "SF_ridge"
                    , "R2_SampCalc"
                    , "R2_mle"
                    , "R2_uniform"
                    , "R2_bootstrap"
                    , "R2_Firths"
                    , "R2_lasso"
                    , "R2_ridge"
                    , "CITL_mle", "CITLLower_mle", "CITLUpper_mle"
                    , "CITL_uniform", "CITLLower_uniform", "CITLUpper_uniform"
                    , "CITL_bootstrap", "CITLLower_bootstrap", "CITLUpper_bootstrap"
                    , "CITL_Firths", "CITLLower_Firths", "CITLUpper_Firths"
                    , "CITL_lasso", "CITLLower_lasso", "CITLUpper_lasso"
                    , "CITL_ridge", "CITLLower_ridge", "CITLUpper_ridge"
                    , "CalSlope_mle", "CalSlopeLower_mle", "CalSlopeUpper_mle"
                    , "CalSlope_uniform", "CalSlopeLower_uniform", "CalSlopeUpper_uniform"
                    , "CalSlope_bootstrap", "CalSlopeLower_bootstrap", "CalSlopeUpper_bootstrap"
                    , "CalSlope_Firths", "CalSlopeLower_Firths", "CalSlopeUpper_Firths"
                    , "CalSlope_lasso", "CalSlopeLower_lasso", "CalSlopeUpper_lasso"
                    , "CalSlope_ridge", "CalSlopeLower_ridge", "CalSlopeUpper_ridge"
                    , "AUC_mle", "AUCLower_mle", "AUCUpper_mle"
                    , "AUC_uniform", "AUCLower_uniform", "AUCUpper_uniform"
                    , "AUC_bootstrap", "AUCLower_bootstrap", "AUCUpper_bootstrap"
                    , "AUC_Firths", "AUCLower_Firths", "AUCUpper_Firths"
                    , "AUC_lasso", "AUCLower_lasso", "AUCUpper_lasso"
                    , "AUC_ridge", "AUCLower_ridge", "AUCUpper_ridge"
                    )
  results <- matrix(nrow = n_iter, ncol = length(result_names))
  colnames(results) <- result_names
  
  ## Repeat the simulation across n_iter number of iterations
  for (i in 1:n_iter) {
    simulation_results <- simulation_singlerun_fnc(P = P,
                                                   Sigma.mat = Sigma.mat,
                                                   Y_prev = Y_prev,
                                                   beta_true = beta_true,
                                                   R2_based_on_maxR2 = R2_based_on_maxR2)
    
    ## Pull out the simulation results from the current iteration
    results[i,"SampleSize"] <- simulation_results$SampleSize 
    results[i,"SampleSizeCriteria"] <- simulation_results$SampleSizeCriteria
    results[i,"SF_uniform"] <- simulation_results$SF_uniform
    results[i,"SF_bootstrap"] <- simulation_results$SF_bootstrap
    results[i,"SF_lasso"] <- simulation_results$SF_lasso 
    results[i,"SF_ridge"] <- simulation_results$SF_ridge
    results[i,"R2_SampCalc"] <- simulation_results$R2_SampCalc 
    results[i,"R2_mle"] <- simulation_results$R2_mle 
    results[i,"R2_uniform"] <- simulation_results$R2_uniform 
    results[i,"R2_bootstrap"] <- simulation_results$R2_bootstrap 
    results[i,"R2_Firths"] <- simulation_results$R2_Firths 
    results[i,"R2_lasso"] <- simulation_results$R2_lasso 
    results[i,"R2_ridge"] <- simulation_results$R2_ridge 
    results[i,"CITL_mle"] <- simulation_results$CITL_mle 
    results[i,"CITLLower_mle"] <- simulation_results$CITLLower_mle 
    results[i,"CITLUpper_mle"] <- simulation_results$CITLUpper_mle 
    results[i,"CITL_uniform"] <- simulation_results$CITL_uniform 
    results[i,"CITLLower_uniform"] <- simulation_results$CITLLower_uniform 
    results[i,"CITLUpper_uniform"] <- simulation_results$CITLUpper_uniform 
    results[i,"CITL_bootstrap"] <- simulation_results$CITL_bootstrap 
    results[i,"CITLLower_bootstrap"] <- simulation_results$CITLLower_bootstrap 
    results[i,"CITLUpper_bootstrap"] <- simulation_results$CITLUpper_bootstrap
    results[i,"CITL_Firths"] <- simulation_results$CITL_Firths 
    results[i,"CITLLower_Firths"] <- simulation_results$CITLLower_Firths 
    results[i,"CITLUpper_Firths"] <- simulation_results$CITLUpper_Firths
    results[i,"CITL_lasso"] <- simulation_results$CITL_lasso
    results[i,"CITLLower_lasso"] <- simulation_results$CITLLower_lasso 
    results[i,"CITLUpper_lasso"] <- simulation_results$CITLUpper_lasso 
    results[i,"CITL_ridge"] <- simulation_results$CITL_ridge 
    results[i,"CITLLower_ridge"] <- simulation_results$CITLLower_ridge 
    results[i,"CITLUpper_ridge"] <- simulation_results$CITLUpper_ridge 
    results[i,"CalSlope_mle"] <- simulation_results$CalSlope_mle
    results[i,"CalSlopeLower_mle"] <- simulation_results$CalSlopeLower_mle 
    results[i,"CalSlopeUpper_mle"] <- simulation_results$CalSlopeUpper_mle 
    results[i,"CalSlope_uniform"] <- simulation_results$CalSlope_uniform 
    results[i,"CalSlopeLower_uniform"] <- simulation_results$CalSlopeLower_uniform 
    results[i,"CalSlopeUpper_uniform"] <- simulation_results$CalSlopeUpper_uniform
    results[i,"CalSlope_bootstrap"] <- simulation_results$CalSlope_bootstrap
    results[i,"CalSlopeLower_bootstrap"] <- simulation_results$CalSlopeLower_bootstrap 
    results[i,"CalSlopeUpper_bootstrap"] <- simulation_results$CalSlopeUpper_bootstrap 
    results[i,"CalSlope_Firths"] <- simulation_results$CalSlope_Firths 
    results[i,"CalSlopeLower_Firths"] <- simulation_results$CalSlopeLower_Firths 
    results[i,"CalSlopeUpper_Firths"] <- simulation_results$CalSlopeUpper_Firths 
    results[i,"CalSlope_lasso"] <- simulation_results$CalSlope_lasso 
    results[i,"CalSlopeLower_lasso"] <- simulation_results$CalSlopeLower_lasso 
    results[i,"CalSlopeUpper_lasso"] <- simulation_results$CalSlopeUpper_lasso 
    results[i,"CalSlope_ridge"] <- simulation_results$CalSlope_ridge 
    results[i,"CalSlopeLower_ridge"] <- simulation_results$CalSlopeLower_ridge 
    results[i,"CalSlopeUpper_ridge"] <- simulation_results$CalSlopeUpper_ridge 
    results[i,"AUC_mle"] <- simulation_results$AUC_mle
    results[i,"AUCLower_mle"] <- simulation_results$AUCLower_mle 
    results[i,"AUCUpper_mle"] <- simulation_results$AUCUpper_mle 
    results[i,"AUC_uniform"] <- simulation_results$AUC_uniform 
    results[i,"AUCLower_uniform"] <- simulation_results$AUCLower_uniform 
    results[i,"AUCUpper_uniform"] <- simulation_results$AUCUpper_uniform
    results[i,"AUC_bootstrap"] <- simulation_results$AUC_bootstrap 
    results[i,"AUCLower_bootstrap"] <- simulation_results$AUCLower_bootstrap
    results[i,"AUCUpper_bootstrap"] <- simulation_results$AUCUpper_bootstrap
    results[i,"AUC_Firths"] <- simulation_results$AUC_Firths
    results[i,"AUCLower_Firths"] <- simulation_results$AUCLower_Firths 
    results[i,"AUCUpper_Firths"] <- simulation_results$AUCUpper_Firths 
    results[i,"AUC_lasso"] <- simulation_results$AUC_lasso 
    results[i,"AUCLower_lasso"] <- simulation_results$AUCLower_lasso 
    results[i,"AUCUpper_lasso"] <- simulation_results$AUCUpper_lasso 
    results[i,"AUC_ridge"] <- simulation_results$AUC_ridge 
    results[i,"AUCLower_ridge"] <- simulation_results$AUCLower_ridge 
    results[i,"AUCUpper_ridge"] <- simulation_results$AUCUpper_ridge 
  }
  ## Return results from the simulation across all iterations
  return(results)
}
