# Developing Clinical Prediction Models: the importance of penalization methods and quantifying bootstrap variability even when adhering to minimum sample size recommendations
Repo for the paper entitled "Developing clinical prediction models when adhering to minimum sample size recommendations: The importance of quantifying bootstrap variability in tuning parameters and predictive performance", which is published in Statististical Methods in Medical Research (2021, DOI: 10.1177/09622802211046388). 

This paper explores the characteristics of CPM performance metrics, upon validation, of models developed using a range of penalisation methods compared with unpenalised maximum likelihood estimation, in derivation data that satisfy formal sample size criteria (Riley et al.). Additionally, it explores the between-sample variability of performance metrics, upon validation, for models developed using each penalisation approach. 

The repo contains the coding scripts and results from the simulation study described in the paper as follows:
## Code sub-folder
This contains the R scripts used in the simulation study. Additionally, this folder also contains the SQL script used to extract the cohort from the MIMIC-III database, which was used as part of the empirical study described in the paper. Much of the simulation code was run on the computational shared facility (CSF) at the University of Manchester.

## Data sub-folder
This contains a .rmd file of the results of the simulation, as described in the paper. The Data sub-folder also contains the results of the MIMIC-III real-world example, as described in the paper.
