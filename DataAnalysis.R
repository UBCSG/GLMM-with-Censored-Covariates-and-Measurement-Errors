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
setwd(here::here("Chapter3"))

## Create results table 
results_new <- data.frame(matrix(ncol = 5, nrow = 9))
results_naive <- data.frame(matrix(ncol = 5, nrow = 2))
results_empirical <- data.frame(matrix(ncol = 5, nrow = 5))
colnames(results_new) <- colnames(results_naive) <- colnames(results_empirical) <- c("Parameter", "Estimates", "Standard Error", "$z$-value", "$p$-value")
results_new$Parameter <- c("$\\alpha_1$","$\\alpha_2$","$\\alpha_3$", "$\\alpha_4$","$\\alpha_5$", "$\\alpha_6$", "$\\alpha_7$","$\\beta_0$","$\\beta_w$")
results_naive$Parameter <- c("$\\beta_0$","$\\beta_w$")
results_empirical$Parameter <- c("$\\alpha_0$","$\\alpha_1$","$\\alpha_2$","$\\beta_0$","$\\beta_w$")

############ Our joint model ==========
# Read in data for Monolix
data_new = list(dataFile = paste0('Data/Cleaned_Data/full2dat_cleaned.csv'),
                      headerTypes =c("id","time","ignore","cens","limit","regressor","regressor","obsid","observation","regressor"),
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
setIndividualParameterVariability(p1 = TRUE, b1 = TRUE, p2 = TRUE, beta1 = TRUE,beta2 = TRUE,beta3 = TRUE,beta4 = TRUE, gamma0 = TRUE, gamma1 = FALSE)
setIndividualParameterDistribution(p1="normal",b1="normal",p2="normal",beta1="normal",beta2="normal",beta3="normal",beta4="normal", gamma0="normal", gamma1="normal")
setPopulationParameterInformation(p1_pop = list(initialValue = 11), 
                                  b1_pop = list(initialValue = 0.2), 
                                  p2_pop = list(initialValue = 6), 
                                  beta1_pop = list(initialValue = 1.3), 
                                  beta2_pop = list(initialValue = 7.6), 
                                  beta3_pop = list(initialValue = 0.1), 
                                  beta4_pop = list(initialValue = 2.4), 
                                  gamma0_pop = list(initialValue = 8), 
                                  gamma1_pop  = list(initialValue = -1.7))
# run the estimation
setScenario(scenario)
runScenario()

# store the estimates
results_new$Estimates <- getEstimatedPopulationParameters()[1:9]
results_new$`Standard Error` <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:9])
results_new$`$z$-value` <- results_new$Estimates/results_new$`Standard Error`
results_new$`$p$-value` <- 2 * pnorm(abs(results_new$`$z$-value`), lower.tail = FALSE)

# save the results in latex file
write.csv(results_new, "Results/results_new.csv")
latex_table<-xtable(results_new, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,3,3,3,3,3)
print(latex_table, file = "Results/results_new.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,9))

# get individual parameters for plots later
indest <- getEstimatedIndividualParameters()$conditionalMean
write.csv(indest, "Data/Other/indest.csv")

# get correlation of estimates
cor <- getEstimatedPopulationParameters()[18:23]
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
print(latex_table, file = "Results/results_cor.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,1))


############ Naive model ===========
data <- read.csv('Data/Cleaned_Data/full2dat_cleaned.csv')
data_naive <- pivot_wider(data, names_from = name, values_from = observation)
data_naive <- data_naive %>% 
  filter(!is.na(CD4)) %>% 
  mutate(censor = ifelse(is.na(rna), 1, censor)) %>% 
  mutate(rna = ifelse(is.na(rna), 2, rna)) %>% 
  mutate(rna = ifelse(censor == 1, log10(50), rna))

model_naive <- glmmPQL(CD4 ~ rna, random = ~ 1|patid, family=binomial, data=data_naive)
results_naive$Estimates <- summary(model_naive)$coefficients$fixed
results_naive$`Standard Error` <- summary(model_naive)$tTable[,2]
results_naive$`$z$-value` <- results_naive$Estimates/results_naive$`Standard Error`
results_naive$`$p$-value` <- 2 * pnorm(abs(results_naive$`$z$-value`), lower.tail = FALSE)


# save the results in latex file
latex_table<-xtable(results_naive, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,3,3,3,3,3)
print(latex_table, file = "Results/results_naive.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,2))


############ Empirical model ============
data_empirical <- data %>% 
  mutate(observation = ifelse(censor == 1 & name=="rna", log10(50), observation))
write.csv(data_empirical, "Data/Other/data_empirical.csv", row.names = FALSE)

# Read in data for Monolix
data_empirical = list(dataFile = paste0('Data/Other/data_empirical.csv'),
                      headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","obsid","observation","ignore"),
                      observationTypes =list(CD4 = "categorical", rna = "continuous"))
modelFile = paste0('Model/model_empirical.txt')
# create a new project by setting a data set and a structural model
newProject(data = data_empirical, modelFile = modelFile)
# set error model and observation distribution
setErrorModel(list(ylog10rna = "constant"))
setObservationDistribution(ylog10rna= "normal")
# set tasks in scenario
scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, gamma0 = TRUE, gamma1 = FALSE)
setPopulationParameterEstimationSettings(nbexploratoryiterations=1000)
setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal", gamma0="normal", gamma1="normal")
setPopulationParameterInformation(gamma0_pop = list(initialValue = 8), 
                                  gamma1_pop  = list(initialValue = -2.5))
# run the estimation
setScenario(scenario)
runScenario()
# store the estimates
results_empirical$Estimates <- getEstimatedPopulationParameters()[1:5]
results_empirical$`Standard Error` <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:5])
results_empirical$`$z$-value` <- results_empirical$Estimates/results_empirical$`Standard Error`
results_empirical$`$p$-value` <- 2 * pnorm(abs(results_empirical$`$z$-value`), lower.tail = FALSE)

# save the results in latex file
latex_table<-xtable(results_empirical, type = "latex",align=c("cccccc"))
digits(latex_table)<-c(0,3,3,3,3,3)
print(latex_table, file = "Results/results_empirical.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,5))

# get individual parameters for plots later
indest <- getEstimatedIndividualParameters()$saem
write.csv(indest, "Data/Other/indest_empirical.csv")

