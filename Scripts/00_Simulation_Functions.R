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
## Define a function to repeat the simulation across all iterations
####-----------------------------------------------------------------------------------------
simulation_nruns_fnc <- function(n_iter,
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
  
  ## Define an empty variable, which will be used to store the results across all iterations
  results <- NULL
  
  ## Repeat the simulation across n_iter number of iterations
  for (iter in 1:n_iter) {
    simulation_results <- simulation_singlerun_fnc(P = P,
                                                   Sigma.mat = Sigma.mat,
                                                   Y_prev = Y_prev,
                                                   beta_true = beta_true,
                                                   R2_based_on_maxR2 = R2_based_on_maxR2)
    #Save results across iterations:
    results <- results %>%
      bind_rows(as_tibble(simulation_results) %>% 
                  mutate("Iteration" = iter, .before = "Model"))
    
    rm(simulation_results)
    
  }
  ## Return results from the simulation across all iterations
  return(results)
}

####-----------------------------------------------------------------------------------------
## Define a function to run the processes within a single iteration
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
    MaxR2 <- 1-(((Y_prev^(Y_prev))*((1-Y_prev)^(1-Y_prev)))^2)
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
  SampleSizeCriteria <- SampleSize$Max_Criteria
  
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
    
    Raw_LP <- predict(BootMod, newdata = DevelopmentData, type = "link")
    Bootstrap_shrinkage[boot] <- as.numeric(coef(glm(DevelopmentData$Y ~ Raw_LP, 
                                                     family = binomial(link = "logit")))[2])
    
  }
  Bootstrap_shrinkage <- mean(Bootstrap_shrinkage)
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
  #re-estimate the intercept ready for prediction:
  Firths_LP <- as.numeric((DevelopmentData %>%
                             select(starts_with("V")) %>%
                             data.matrix()) %*% (as.numeric(coef(Firths_CPM))[-1]))
  Firths.beta <- c(coef(glm(DevelopmentData$Y ~ offset(Firths_LP), 
                            family = binomial(link = "logit")))[1],
                   (as.numeric(coef(Firths_CPM))[-1])) 
  #Apply shrunk model to the validation set:
  Firths_CPM_validation_LP <- as.numeric(cbind(1,
                                               ValidationData %>%
                                                 select(starts_with("V")) %>%
                                                 data.matrix()) %*% Firths.beta)
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
  PredictivePerformance <- map_dfr(list("MLE" = ValidationData$MLE_CPM_Pi
                                        , "Uniform" = ValidationData$UniformShrunk_CPM_Pi
                                        , "Bootstrap" = ValidationData$BootstrapShrunk_CPM_Pi
                                        , "Firth" = ValidationData$Firths_CPM_Pi
                                        , "LASSO" = ValidationData$LASSO_CPM_Pi
                                        , "Ridge" = ValidationData$Ridge_CPM_Pi),
                                   Performance_fnc,
                                   Response = ValidationData$Y,
                                   .id = "Model") %>%
    #Extract relevant information to store per iteration:
    mutate("SampleSize" = N_Dev,
           "SampleSizeCriteria" = SampleSizeCriteria,
           
           "SF_uniform" = uniform_shrinkage,
           "SF_bootstrap" = Bootstrap_shrinkage,
           "SF_lasso" = LASSO_CPM$lambda.min,
           "SF_ridge" = Ridge_CPM$lambda.min,
           
           "R2_SampCalc" = R2_for_SampCalc,
           .after = "Model")
  
  
  ### Return the results
  return(PredictivePerformance)
}


####-----------------------------------------------------------------------------------------
## Define function to assess predictive performance of each model
####-----------------------------------------------------------------------------------------
Performance_fnc <- function(PredictedRisk, Response) {
  LP <- log(PredictedRisk / (1 - PredictedRisk))
  
  CITL_mod <- glm(Response ~ offset(LP), family = binomial(link = "logit"))
  CITL <- as.numeric(coef(CITL_mod)[1])
  CITLSE <- sqrt(vcov(CITL_mod)[1,1])  
  
  CalSlope_mod <- glm(Response ~ LP, family = binomial(link = "logit"))
  CalSlope <- as.numeric(coef(CalSlope_mod)[2])
  CalSlopeSE <- sqrt(vcov(CalSlope_mod)[2,2])
  
  AUC <- roc(response = Response,
             predictor = PredictedRisk,
             direction = "<",
             levels = c(0,1),
             ci = TRUE)
  
  LR <- -2 * (as.numeric(logLik(glm(Response ~ 1, family = binomial(link = "logit")))) - 
                as.numeric(logLik(CITL_mod)))
  R2_coxsnell <- 1 - exp(-LR/length(Response))
  
  return(list("CITL" = CITL,
              "CITLSE" = CITLSE,
              "CalSlope" = CalSlope,
              "CalSlopeSE" = CalSlopeSE,
              "AUC" = as.numeric(AUC$auc),
              "AUCSE" = sqrt(var(AUC)),
              "R2" = R2_coxsnell))
}

