library(lixoftConnectors)
library(mvtnorm)
library(nlme)
library(MASS)
library(saemix)
library(tidyverse)
library(testit)
library(here)
library(docopt)
library(ggplot2)
library(lmtest)
library(gridExtra)
library(splines)
library(stats)
library(GLDEX)
library(lme4)
library(xtable)
library(berryFunctions)
initializeLixoftConnectors(software="monolix")

## Create results table 
results_new <- data.frame(matrix(ncol = 5, nrow = 19))
results_naive <- data.frame(matrix(ncol = 5, nrow = 3))
results_empirical <- data.frame(matrix(ncol = 5, nrow = 11))
colnames(results_new) <- colnames(results_naive) <- colnames(results_empirical) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
results_new$Parameter <- c("$\\alpha_1$","$\\alpha_2$","$\\alpha_{21}$", "$\\alpha_3$", "$\\alpha_4$","$\\alpha_5$", "$\\alpha_6$", "$\\alpha_7$","$\\beta_0$","$\\beta_w$","$A_{11}$","$A_{22}$","$A_{33}$","$A_{44}$","$A_{55}$","$A_{66}$","$A_{77}$","$\\sigma_1$","$b$")
results_naive$Parameter <- c("$\\beta_0$","$\\beta_w$","$b$")
results_empirical$Parameter <- c("$\\alpha_0$","\\gamma","$\\alpha_1$","$\\alpha_2$","$\\beta_0$","$\\beta_w$","$D_{11}$","$D_{22}$","$D_{33}$", "$\\sigma_2$", "$b$")

############ Our joint model ==========
# Read in data for Monolix
data_new = list(dataFile = paste0('Data/Cleaned_Data/full2dat_cleaned_with_cov.csv'),
                      headerTypes =c("id","time","ignore","cens","limit","regressor","regressor","ignore", "ignore","ignore", "contcov","contcov","contcov","contcov", "obsid","observation","regressor"),
                      observationTypes =list(CD4 = "categorical", rna = "continuous"))
modelFile = paste0('Model/model_joint.txt')
# create a new project by setting a data set and a structural model
newProject(data = data_new, modelFile = modelFile)
# set error model and observation distribution
setErrorModel(list(yrna = "constant"))
setObservationDistribution(yrna= "normal")
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setCorrelationBlocks(ID = list(c("b1", "beta1", "beta2", "beta3")))
setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE, beta1 = TRUE, beta2 = TRUE,beta3 = TRUE,beta4 = TRUE, gamma0 = TRUE, gamma1 = FALSE)
setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",beta1="normal",beta2="normal",beta3="normal",beta4="normal", gamma0="normal", gamma1="normal")
setPopulationParameterInformation(p1_pop = list(initialValue = 11), 
                                  b1_pop = list(initialValue = 0.3), 
                                  p2_pop = list(initialValue = 5.7), 
                                  beta1_pop = list(initialValue = 1.4), 
                                  beta2_pop = list(initialValue = 5), 
                                  beta3_pop = list(initialValue = 0.2), 
                                  beta4_pop = list(initialValue = 2.4), 
                                  gamma0_pop = list(initialValue = 9.7), 
                                  gamma1_pop  = list(initialValue = -2))
setCovariateModel(b1 = c(age_st = TRUE))
# run the estimation
setScenario(scenario)
runScenario()

# store the estimates
results_new$Estimates[1:10] <- getEstimatedPopulationParameters()[1:10]
results_new$Estimates[11:17] <- (getEstimatedPopulationParameters()[11:17])^2
results_new$`Standard Error`[1:10] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:10])
results_new$`$z$-value`[1:10] <- results_new$Estimates[1:10]/results_new$`Standard Error`[1:10]
results_new$`$p$-value`[1:10] <- 2 * pnorm(abs(results_new$`$z$-value`[1:10]), lower.tail = FALSE)
results_new$Estimates[18] <- getEstimatedPopulationParameters()[25]
results_new$Estimates[19] <- getEstimatedPopulationParameters()[18]

# save the results in latex file
write.csv(results_new, "Results/results_with_cov_new.csv")
latex_table<-xtable(results_new, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,3,3,3,3,3)
print(latex_table, file = "Results/results_with_cov_new.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,19))

# get individual parameters for plots later
indest <- getEstimatedIndividualParameters()$conditionalMean
write.csv(indest, "Data/Other/indest_with_cov.csv")

# get correlation of estimates
cor <- getEstimatedPopulationParameters()[19:24]
cortable <- as.data.frame(matrix(NA, 4, 5))
colnames(cortable) <- c(" ", "\\alpha_2", "\\alpha_4", "\\alpha_5", "\\alpha_6")
cortable$` ` <- c("\\alpha_2", "\\alpha_4", "\\alpha_5", "\\alpha_6")
cortable$`\\alpha_2` <- c(1, cor[1:3])
cortable$`\\alpha_4` <- c(cor[1], 1, cor[4:5])
cortable$`\\alpha_5`<- c(cor[2], cor[4], 1, cor[6])
cortable$`\\alpha_6`<- c(cor[3], cor[5:6], 1)
# save the results in latex file
latex_table <- xtable(cortable, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,0,3,3, 3,3)
print(latex_table, file = "Results/results_with_cov_cor.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,1))


############ Naive model ===========
data <- read.csv('Data/Cleaned_Data/full2dat_cleaned_with_cov.csv')
data_naive <- pivot_wider(data, names_from = name, values_from = observation)
data_naive <- data_naive %>% 
  filter(!is.na(CD4)) %>% 
  mutate(censor = ifelse(is.na(rna), 1, censor)) %>% 
  mutate(rna = ifelse(is.na(rna), 2, rna)) %>% 
  mutate(rna = ifelse(censor == 1, log10(50), rna))
write.csv(data_naive, "Data/Other/data_naive_with_cov.csv", row.names = FALSE)

# Read in data for Monolix
data_naive = list(dataFile = paste0('Data/Other/data_naive_with_cov.csv'),
                  headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","ignore","ignore","ignore","ignore","ignore","ignore","ignore","ignore","contcov","observation"),
                  observationTypes =list(CD4 = "categorical"))
modelFile = paste0('Model/model_naive.txt')
# create a new project by setting a data set and a structural model
newProject(data = data_naive, modelFile = modelFile)
setCovariateModel(gamma0 = c(rna=TRUE))
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setIndividualParameterVariability(gamma0 = TRUE)
setIndividualParameterDistribution(gamma0="normal")
setPopulationParameterInformation(gamma0_pop = list(initialValue = 9.7))
# run the estimation
setScenario(scenario)
runScenario()
# store the estimates
results_naive$Estimates <- getEstimatedPopulationParameters()[1:3]
results_naive$`Standard Error`[1:2] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:2])
results_naive$`$z$-value`[1:2] <- results_naive$Estimates[1:2]/results_naive$`Standard Error`[1:2]
results_naive$`$p$-value`[1:2] <- 2 * pnorm(abs(results_naive$`$z$-value`[1:2]), lower.tail = FALSE)

# save the results in latex file
latex_table<-xtable(results_naive, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,3,3,3,3,3)
print(latex_table, file = "Results/results_with_cov_naive.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,3))


############ Empirical model ============
data_empirical <- data %>% 
  mutate(observation = ifelse(censor == 1 & name=="rna", log10(50), observation)) %>% 
  mutate(day = day/30)
write.csv(data_empirical, "Data/Other/data_empirical.csv", row.names = FALSE)

# Read in data for Monolix
data_empirical = list(dataFile = paste0('Data/Other/data_empirical.csv'),
                      headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","ignore", "ignore","ignore", "contcov","contcov","contcov","contcov","obsid","observation","ignore"),
                      observationTypes =list(CD4 = "categorical", rna = "continuous"))
modelFile = paste0('Model/model_empirical.txt')
# create a new project by setting a data set and a structural model
newProject(data = data_empirical, modelFile = modelFile)
# set error model and observation distribution
setErrorModel(list(yrna = "constant"))
setObservationDistribution(yrna= "normal")
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, gamma0 = TRUE, gamma1 = FALSE)
setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal", gamma0="normal", gamma1="normal")
setPopulationParameterInformation(gamma0_pop = list(initialValue = 9.7), 
                                  gamma1_pop  = list(initialValue = -2))
setCovariateModel(alpha0 = c(age_st = TRUE))
# run the estimation
setScenario(scenario)
runScenario()
# store the estimates
results_empirical$Estimates[1:6] <- getEstimatedPopulationParameters()[1:6]
results_empirical$Estimates[7:9] <- (getEstimatedPopulationParameters()[7:9])^2
results_empirical$`Standard Error`[1:6] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:6])
results_empirical$`$z$-value`[1:6] <- results_empirical$Estimates[1:6]/results_empirical$`Standard Error`[1:6]
results_empirical$`$p$-value`[1:6] <- 2 * pnorm(abs(results_empirical$`$z$-value`[1:6]), lower.tail = FALSE)
results_empirical$Estimates[10] <- getEstimatedPopulationParameters()[11]
results_empirical$Estimates[11] <- getEstimatedPopulationParameters()[10]

# save the results in latex file
write.csv(results_empirical, "Results/results_with_cov_empirical.csv")
latex_table<-xtable(results_empirical, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,3,3,3,3,3)
print(latex_table, file = "Results/results_with_cov_empirical.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,11))

# get individual parameters for plots later
indest <- getEstimatedIndividualParameters()$saem
write.csv(indest, "Data/Other/indest_with_cov_empirical.csv")

