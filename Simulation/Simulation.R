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
setwd(here::here())

num_sim <- 100  # number of simulation
num_patient <- 100  # number of patients
censor_value <- log10(25)  # censor value 
num_boots = 200 # number of bootstrap iterations


# Set true values for alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, sigma, 
# and covariance matrix A for viral loads before ART interruption and following ART interruption
# Follow model: w_ij=[log10(exp(alpha_{1i}-alpha_{2i}*t_{ij})+exp(alpha_3i))]I(t_{ij} < T_i) + 
#         [alpha_{4i}t_{ik}/(t_{ik}+exp(alpha_{5i}-alpha_{6i}t_{ik}))+alpha_{7i}]I(t_{ik} >= T_i)+e_{ij}
# w_{ij} = w_{ij}^* + e_{ij}
# Random effects: alpha_{il} = alpha_l + a_{li}, l = 1,..,7
# where e_{ij}~N(0,sigma^2) and a_i~N(0,A)
alpha1 = 11
alpha2 = 0.3
alpha3 = 2.8
alpha4 = 1.8
alpha5 = 8
alpha6 = 0.4
alpha7 = 1.3
sigma <- 0.5
A <- diag(c(2, 0.01, 0.05, 0.05, 0, 0.01, 0))


# Set true value for beta0, beta1, and B for binary CD4 model
# Follow model: log(P(y_{ij}=1)/P(y_{ij}=0))=(beta0 + b_{0i}) + beta1 * w_{ij}^*
beta0 <- 8
beta1 <- -2.5
B <- 0.1

parnames <- c("beta0","beta1")
parnames_viral <- c("alpha1", "alpha2", "alpha3", "alpha4","alpha5", "alpha6","alpha7")

# Choose time from true_value dataset. 
time0 <- c(0, 8, 14, 28, 56, 65, 84, 87, 90 , 93, 96,  99, 102, 104, 108, 111, 116, 120, 125,140) 


#Create a data frame for simulation estimates.
ncol <- length(parnames)
ncol_viral <- length(parnames_viral)

estimates_naive <- estimates_empirical <- estimates_new <- data.frame(matrix(nrow=num_sim,ncol=ncol))
colnames(estimates_naive) <- colnames(estimates_empirical) <- colnames(estimates_new) <- parnames
estimates_new_viral <- data.frame(matrix(nrow=num_sim,ncol=ncol_viral))
colnames(estimates_new_viral) <- parnames_viral

SE_new <- SE_empirical <- SE_naive <- data.frame(matrix(nrow=num_sim,ncol=ncol))
colnames(SE_naive) <- colnames(SE_empirical) <- colnames(SE_new) <- parnames
SE_new_viral <- data.frame(matrix(nrow=num_sim,ncol=ncol_viral))
colnames(SE_new_viral) <- parnames_viral


bootstrap_nonpar_new<- bootstrap_par_new <- bootstrap_par_naive <- bootstrap_par_empirical <- data.frame(matrix(nrow=num_boots,ncol=ncol))
colnames(bootstrap_par_new) <- colnames(bootstrap_nonpar_new) <- colnames(bootstrap_par_naive) <- colnames(bootstrap_par_empirical) <- parnames
bootstrap_par_new_viral <- bootstrap_nonpar_new_viral <- data.frame(matrix(nrow=num_boots,ncol =ncol_viral))
colnames(bootstrap_par_new_viral) <- colnames(bootstrap_par_new_viral) <- parnames_viral


bootstrap_par_SE_new <- bootstrap_par_SE_naive <- bootstrap_par_SE_empirical <- data.frame(matrix(nrow = num_sim,ncol = ncol))
colnames(bootstrap_par_SE_new) <- colnames(bootstrap_par_SE_naive) <- colnames(bootstrap_par_SE_empirical) <- parnames
bootstrap_par_SE_new_viral<- data.frame(matrix(nrow=num_sim,ncol = ncol_viral))
colnames(bootstrap_par_SE_new_viral) <- parnames_viral


true_value <- c(beta0, beta1)
true_value_viral <- c(alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7)

i=1
for (i in 1:100){
  print(glue::glue("i = ",i))
  
  # Create dataset for viral load
  PATIENT <- rep(1:num_patient, each = length(time0))
  time <- rep(time0, num_patient)
  decay <- rep(c(rep(1, 7), rep(0, 13)),num_patient)
  rebound <- rep(c(rep(0, 7), rep(1, 13)),num_patient)
  Ti <- rep(84, num_patient*20)
  data0 <- as.data.frame(cbind(PATIENT, time, decay, rebound, Ti))
  
  # simulate a_i, b_i, and e_ij
  ai_sim <- mvrnorm(num_patient, rep(0,dim(A)[1]), A)
  colnames(ai_sim) <- c("a1_sim","a2_sim","a3_sim","a4_sim","a5_sim","a6_sim","a7_sim")
  b0_sim <- rnorm(num_patient, B)
  ranef_sim <- cbind(PATIENT=c(1:num_patient), ai_sim, b0_sim)
  data <- merge(data0, ranef_sim, by="PATIENT", all = TRUE)
  data <- data %>% 
    mutate(e = rnorm(nrow(data), 0, sigma)) %>% 
    mutate(w_star = log10(exp(alpha1+a1_sim-(alpha2+a2_sim)*time)+exp(alpha3+a3_sim))*decay + ((alpha4+a4_sim) * (time-Ti) / ((time-Ti)  + exp(alpha5 + a5_sim - (alpha6+a6_sim) * (time-Ti))) + (alpha7 + a7_sim))*rebound) %>% 
    mutate(log10rna = w_star + e) %>% 
    mutate(censor = ifelse(log10rna < censor_value, 1, 0)) %>% 
    mutate(log10rna = ifelse(log10rna < censor_value, censor_value, log10rna)) %>% 
    mutate(limit = 0) %>% 
    mutate(xb = beta0 + b0_sim + beta1 * w_star) %>% 
    mutate(p = 1/(1+exp(-xb))) %>% 
    mutate(CD4 = rbinom(n = nrow(data), size = 1, p = p)) 
  data <- data %>% 
    dplyr::select(c(PATIENT, time, decay, rebound, Ti, censor, limit, log10rna, CD4))
  data_new <- pivot_longer(data, cols = c(8:9), names_to = "name", values_to = "observation")
  data_new <- data_new %>% 
    mutate(censor = ifelse(name == "CD4", 0, censor))
  data_empirical <- data_new %>% 
    mutate(observation = ifelse(censor == 1 & name=="log10rna", 0.5*observation, observation))
  data_naive <- data %>% 
    mutate(log10rna = ifelse(censor == 1, 0.5*log10rna, log10rna))
  summary(as.factor(data$censor))
  summary(as.factor(data$CD4))
  ggplot(data[PATIENT%in%c(1:5),], aes(x=time, y=log10rna)) +
    geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
    geom_line(aes(group=PATIENT)) +
    scale_x_continuous("Day") +
    scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
    scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
    scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
    labs(color = "ART status", fill = "Data type")+ggtitle("Plot for all observations")+
    theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
  
  # Save the simulated data locally
  write.table(data_new, "Simulation/SimData/sim_data_new.txt", sep = ",", quote = FALSE, row.names = FALSE)
  write.table(data_naive, "Simulation/SimData/sim_data_naive.txt", sep = ",", quote = FALSE, row.names = FALSE)
  write.table(data_empirical, "Simulation/SimData/sim_data_empirical.txt", sep = ",", quote = FALSE, row.names = FALSE)
  
  # Our model 
  # Read in data for Monolix
  data_joint = list(dataFile = paste0('Simulation/SimData/sim_data_new.txt'),
                    headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                    observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
  modelFile = paste0('Model/model_joint.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data_joint, modelFile = modelFile)

  setErrorModel(ylog10rna = "constant")
  setObservationDistribution(ylog10rna= "normal")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(alpha1 = TRUE, alpha2 = TRUE, alpha3 = TRUE, alpha4 = TRUE, alpha5 = FALSE, alpha6 = TRUE, alpha7 = FALSE, beta0 = TRUE, beta1 = FALSE)
  setIndividualParameterDistribution(alpha1="normal",alpha2="normal",alpha3="normal",alpha4="normal",alpha5="normal",alpha6="normal",alpha7="normal", beta0="normal", beta1="normal")
  setPopulationParameterInformation(alpha1_pop = list(initialValue = 11), 
                                    alpha2_pop = list(initialValue = 0.3),
                                    alpha3_pop = list(initialValue = 2.8),
                                    alpha4_pop = list(initialValue = 1.8),
                                    alpha5_pop = list(initialValue = 8),
                                    alpha6_pop = list(initialValue = 0.4),
                                    alpha7_pop = list(initialValue = 1.3),
                                    beta0_pop = list(initialValue = 8), 
                                    beta1_pop  = list(initialValue = -2.5))
  # run the estimation
  setScenario(scenario)
  runScenario()
  # store the estimates
  estimates_new[i,] <- getEstimatedPopulationParameters()[8:9]
  SE_new[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][8:9])
  estimates_new_viral[i,] <- getEstimatedPopulationParameters()[1:7]
  SE_new_viral[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:7])
  
  # Naive model
  # Read in data for Monolix
  data_naive = list(dataFile = paste0('Simulation/SimData/sim_data_naive.txt'),
                        headerTypes =c("id","ignore","ignore","ignore","ignore","ignore","ignore","contcov","observation"),
                        observationTypes =list(CD4 = "categorical"))
  modelFile = paste0('Model/model_naive.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data_naive, modelFile = modelFile)
  
  # covariates on beta0
  setCovariateModel(beta0 = c(log10rna=TRUE))
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(beta0 = TRUE)
  setIndividualParameterDistribution(beta0="normal")
  setPopulationParameterInformation(beta0_pop = list(initialValue = 8))
  # run the estimation
  setScenario(scenario)
  runScenario()
  # store the estimates
  estimates_naive[i,] <- getEstimatedPopulationParameters()[1:2]
  SE_naive[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:2])
  
  # Empirical model
  # Read in data for Monolix
  data_empirical = list(dataFile = paste0('Simulation/SimData/sim_data_empirical.txt'),
                        headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","obsid","observation"),
                        observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
  modelFile = paste0('Model/model_empirical.txt')
  # create a new project by setting a data set and a structural model
  newProject(data = data_empirical, modelFile = modelFile)
  setErrorModel(list(ylog10rna = "constant"))
  setObservationDistribution(ylog10rna = "normal")
  # set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, beta0 = TRUE, beta1 = FALSE)
  setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal", beta0="normal", beta1="normal")
  setPopulationParameterInformation(beta0_pop = list(initialValue = 8), 
                                    beta1_pop  = list(initialValue = -2.5))
  # run the estimation
  setScenario(scenario)
  runScenario()
  # store the estimates
  estimates_empirical[i,] <- getEstimatedPopulationParameters()[4:5]
  SE_empirical[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][4:5])
  
  while(any(is.na(c(SE_empirical[i, ], SE_new[i, ], SE_naive[i, ], SE_new_viral[i, ])))){
    # simulate a_i, b_i, and e_ij
    ai_sim <- mvrnorm(num_patient, rep(0,dim(A)[1]), A)
    colnames(ai_sim) <- c("a1_sim","a2_sim","a3_sim","a4_sim","a5_sim","a6_sim","a7_sim")
    b0_sim <- rnorm(num_patient, B)
    ranef_sim <- cbind(PATIENT=c(1:num_patient), ai_sim, b0_sim)
    data <- merge(data0, ranef_sim, by="PATIENT", all = TRUE)
    data <- data %>% 
      mutate(e = rnorm(nrow(data), 0, sigma)) %>% 
      mutate(w_star = log10(exp(alpha1+a1_sim-(alpha2+a2_sim)*time)+exp(alpha3+a3_sim))*decay + ((alpha4+a4_sim) * (time-Ti) / ((time-Ti)  + exp(alpha5 + a5_sim - (alpha6+a6_sim) * (time-Ti))) + (alpha7 + a7_sim))*rebound) %>% 
      mutate(log10rna = w_star + e) %>% 
      mutate(censor = ifelse(log10rna < censor_value, 1, 0)) %>% 
      mutate(log10rna = ifelse(log10rna < censor_value, censor_value, log10rna)) %>% 
      mutate(limit = 0) %>% 
      mutate(xb = beta0 + b0_sim + beta1 * w_star) %>% 
      mutate(p = 1/(1+exp(-xb))) %>% 
      mutate(CD4 = rbinom(n = nrow(data), size = 1, p = p)) 
    data <- data %>% 
      dplyr::select(c(PATIENT, time, decay, rebound, Ti, censor, limit, log10rna, CD4))
    data_new <- pivot_longer(data, cols = c(8:9), names_to = "name", values_to = "observation")
    data_new <- data_new %>% 
      mutate(censor = ifelse(name == "CD4", 0, censor))
    data_empirical <- data_new %>% 
      mutate(observation = ifelse(censor == 1 & name=="log10rna", 0.5*observation, observation))
    data_naive <- data %>% 
      mutate(log10rna = ifelse(censor == 1, 0.5*log10rna, log10rna))

    # Save the simulated data locally
    write.table(data_new, "Simulation/SimData/sim_data_new.txt", sep = ",", quote = FALSE, row.names = FALSE)
    write.table(data_naive, "Simulation/SimData/sim_data_naive.txt", sep = ",", quote = FALSE, row.names = FALSE)
    write.table(data_empirical, "Simulation/SimData/sim_data_empirical.txt", sep = ",", quote = FALSE, row.names = FALSE)
    
    # Our model 
    # Read in data for Monolix
    data_joint = list(dataFile = paste0('Simulation/SimData/sim_data_new.txt'),
                      headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                      observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Model/model_joint.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_joint, modelFile = modelFile)
    
    setErrorModel(ylog10rna = "constant")
    setObservationDistribution(ylog10rna= "normal")
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(alpha1 = TRUE, alpha2 = TRUE, alpha3 = TRUE, alpha4 = TRUE, alpha5 = FALSE, alpha6 = TRUE, alpha7 = FALSE, beta0 = TRUE, beta1 = FALSE)
    setIndividualParameterDistribution(alpha1="normal",alpha2="normal",alpha3="normal",alpha4="normal",alpha5="normal",alpha6="normal",alpha7="normal", beta0="normal", beta1="normal")
    setPopulationParameterInformation(alpha1_pop = list(initialValue = 11), 
                                      alpha2_pop = list(initialValue = 0.3),
                                      alpha3_pop = list(initialValue = 2.8),
                                      alpha4_pop = list(initialValue = 1.8),
                                      alpha5_pop = list(initialValue = 8),
                                      alpha6_pop = list(initialValue = 0.4),
                                      alpha7_pop = list(initialValue = 1.3),
                                      beta0_pop = list(initialValue = 8), 
                                      beta1_pop  = list(initialValue = -2.5))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    estimates_new[i,] <- getEstimatedPopulationParameters()[8:9]
    SE_new[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][8:9])
    estimates_new_viral[i,] <- getEstimatedPopulationParameters()[1:7]
    SE_new_viral[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:7])
    
    # Read in data for Monolix
    data_naive = list(dataFile = paste0('Simulation/SimData/sim_data_naive.txt'),
                      headerTypes =c("id","ignore","ignore","ignore","ignore","ignore","ignore","contcov","observation"),
                      observationTypes =list(CD4 = "categorical"))
    modelFile = paste0('Model/model_naive.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_naive, modelFile = modelFile)

    setCovariateModel(beta0 = c(log10rna=TRUE))
    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(beta0 = TRUE)
    setIndividualParameterDistribution(beta0="normal")
    setPopulationParameterInformation(beta0_pop = list(initialValue = 8))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    estimates_naive[i,] <- getEstimatedPopulationParameters()[1:2]
    SE_naive[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][1:2])
    
    
    # Empirical model
    # Read in data for Monolix
    data_empirical = list(dataFile = paste0('Simulation/SimData/sim_data_empirical.txt'),
                          headerTypes =c("id","time","ignore","ignore","ignore","ignore","ignore","obsid","observation"),
                          observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Model/model_empirical.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_empirical, modelFile = modelFile)
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
    setIndividualParameterVariability(alpha0 = TRUE, alpha1 = TRUE, alpha2 = TRUE, beta0 = TRUE, beta1 = FALSE)
    setIndividualParameterDistribution(alpha0="normal",alpha1="normal",alpha2="normal", beta0="normal", beta1="normal")
    setPopulationParameterInformation(beta0_pop = list(initialValue = 8), 
                                      beta1_pop  = list(initialValue = -2.5))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    estimates_empirical[i,] <- getEstimatedPopulationParameters()[4:5]
    SE_empirical[i,] <- as.numeric(getEstimatedStandardErrors()[["stochasticApproximation"]][["se"]][4:5])
  }
  
  for (j in 1:num_boots){
    print(glue::glue("i = ",i," j = ",j))
    
    # Parametric Bootstrap
    # simulate a_i, b_i, and e_ij
    ai_sim <- mvrnorm(num_patient, rep(0,dim(A)[1]), A)
    colnames(ai_sim) <- c("a1_sim","a2_sim","a3_sim","a4_sim","a5_sim","a6_sim","a7_sim")
    b0_sim <- rnorm(num_patient, B)
    ranef_sim <- cbind(PATIENT=c(1:num_patient), ai_sim, b0_sim)
    data <- merge(data0, ranef_sim, by="PATIENT", all = TRUE)
    data <- data %>% 
      mutate(e = rnorm(nrow(data), 0, sigma)) %>% 
      mutate(w_star = log10(exp(estimates_new_viral[i, 1]+a1_sim-(estimates_new_viral[i, 2]+a2_sim)*time)+exp(estimates_new_viral[i, 3]+a3_sim))*decay + ((estimates_new_viral[i, 4]+a4_sim) * (time-Ti) / ((time-Ti)  + exp(estimates_new_viral[i, 5] + a5_sim - (estimates_new_viral[i, 6]+a6_sim) * (time-Ti))) + estimates_new_viral[i, 7] + a7_sim)*rebound) %>% 
      mutate(log10rna = w_star + e) %>% 
      mutate(censor = ifelse(log10rna < censor_value, 1, 0)) %>% 
      mutate(log10rna = ifelse(log10rna < censor_value, censor_value, log10rna)) %>% 
      mutate(limit = 0) %>% 
      mutate(xb = estimates_new[i, 1] + b0_sim + estimates_new[i, 2] * w_star) %>% 
      mutate(p = 1/(1+exp(-xb))) %>% 
      mutate(CD4 = rbinom(n = nrow(data), size = 1, p = p)) 
    data <- data %>% 
      dplyr::select(c(PATIENT, time, decay, rebound, Ti, censor, limit, log10rna, CD4))
    data_new_boot <- pivot_longer(data, cols = c(8:9), names_to = "name", values_to = "observation")
    data_new_boot <- data_new_boot %>% 
      mutate(censor = ifelse(name == "CD4", 0, censor))
    ggplot(data[PATIENT%in%c(1:5),], aes(x=time, y=log10rna)) +
      geom_point(aes(fill=factor(censor)),size=2, shape=21, stroke=0) +
      geom_line(aes(group=PATIENT)) +
      scale_x_continuous("Day") +
      scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"))+
      scale_fill_manual(values=c("black","red"),labels = c("Observed data", "Censored data"))+
      scale_color_manual(values=c("steelblue","black"),labels = c("Following ART interruption", "During ART"))+
      labs(color = "ART status", fill = "Data type")+ggtitle("Plot for all observations")+
      theme_classic()+theme(text=element_text(size=14),panel.grid.major=element_line(colour = "grey90"), panel.grid.minor =element_line(colour = "grey90"))
    
    data_empirical_boot <- data_new_boot %>% 
      mutate(observation = ifelse(censor == 1 & name=="log10rna", 0.5*observation, observation))
    data_naive_boot <- data %>% 
      mutate(log10rna = ifelse(censor == 1, 0.5*log10rna, log10rna))
    
    # Save the simulated data locally
    write.table(data_new_boot, "Simulation/SimData/sim_data_new_boot.txt", sep = ",", quote = FALSE, row.names = FALSE)
    write.table(data_empirical_boot, "Simulation/SimData/sim_data_empirical_boot.txt", sep = ",", quote = FALSE, row.names = FALSE)
    write.table(data_naive_boot, "Simulation/SimData/sim_data_naive_boot.txt", sep = ",", quote = FALSE, row.names = FALSE)
    
    # Our model
    # Read in data for Monolix
    data_new_boot = list(dataFile = paste0('Simulation/SimData/sim_data_new_boot.txt'),
                         headerTypes =c("id","time","regressor", "regressor", "regressor","cens","limit","obsid","observation"),
                         observationTypes =list(CD4 = "categorical", log10rna = "continuous"))
    modelFile = paste0('Model/model_joint.txt')
    # create a new project by setting a data set and a structural model
    newProject(data = data_new_boot, modelFile = modelFile)

    setErrorModel(ylog10rna = "constant")
    setObservationDistribution(ylog10rna= "normal")

    # set tasks in scenario
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setIndividualParameterVariability(alpha1 = TRUE, alpha2 = TRUE, alpha3 = TRUE, alpha4 = TRUE, alpha5 = FALSE, alpha6 = TRUE, alpha7 = FALSE, beta0 = TRUE, beta1 = FALSE)
    setIndividualParameterDistribution(alpha1="normal",alpha2="normal",alpha3="normal",alpha4="normal",alpha5="normal",alpha6="normal",alpha7="normal", beta0="normal", beta1="normal")
    setPopulationParameterInformation(alpha1_pop = list(initialValue = estimates_new_viral[i, 1]),
                                      alpha2_pop = list(initialValue = estimates_new_viral[i, 2]),
                                      alpha3_pop = list(initialValue = estimates_new_viral[i, 3]),
                                      alpha4_pop = list(initialValue = estimates_new_viral[i, 4]),
                                      alpha5_pop = list(initialValue = estimates_new_viral[i, 5]),
                                      alpha6_pop = list(initialValue = estimates_new_viral[i, 6]),
                                      alpha7_pop = list(initialValue = estimates_new_viral[i, 7]),
                                      beta0_pop = list(initialValue = estimates_new[i, 1]),
                                      beta1_pop  = list(initialValue = estimates_new[i, 2]))
    # run the estimation
    setScenario(scenario)
    runScenario()
    # store the estimates
    bootstrap_par_new[j,] <- getEstimatedPopulationParameters()[8:9]
    bootstrap_par_new_viral[j,] <- getEstimatedPopulationParameters()[1:7]
  }
  bootstrap_par_SE_new[i,] <- sapply(bootstrap_par_new,sd)
  bootstrap_par_SE_new_viral[i,] <- sapply(bootstrap_par_new_viral,sd)
  
  saveRDS(estimates_new, "Simulation/SimData/estimates_new_boot.RDS")
  saveRDS(estimates_new_viral, "Simulation/SimData/estimates_new_viral_boot.RDS")
  saveRDS(estimates_naive, "Simulation/SimData/estimates_naive_boot.RDS")
  saveRDS(estimates_empirical, "Simulation/SimData/estimates_empirical_boot.RDS")
  saveRDS(SE_new, "Simulation/SimData/SE_new_boot.RDS")
  saveRDS(SE_new_viral, "Simulation/SimData/SE_new_viral_boot.RDS")
  saveRDS(SE_naive, "Simulation/SimData/SE_naive_boot.RDS")
  saveRDS(SE_empirical, "Simulation/SimData/SE_empirical_boot.RDS")
  saveRDS(bootstrap_par_SE_new,"Simulation/SimData/bootstrap_par_SE_new_boot.RDS")
  saveRDS(bootstrap_par_SE_new_viral,"Simulation/SimData/bootstrap_par_SE_new_viral_boot.RDS")
}  



estimates_new <- readRDS("Simulation/SimData/estimates_new_boot.RDS")
SE_new <- readRDS("Simulation/SimData/SE_new_boot.RDS")
bootstrap_par_SE_new  <- readRDS("Simulation/SimData/bootstrap_par_SE_new_boot.RDS")
bootstrap_par_SE_naive  <- readRDS("Simulation/SimData/bootstrap_par_SE_naive_boot.RDS")
bootstrap_par_SE_empirical  <- readRDS("Simulation/SimData/bootstrap_par_SE_empirical_boot.RDS")
estimates_new_viral <- readRDS("Simulation/SimData/estimates_new_viral_boot.RDS")
SE_new_viral <- readRDS("Simulation/SimData/SE_new_viral_boot.RDS")
bootstrap_par_SE_new_viral  <- readRDS("Simulation/SimData/bootstrap_par_SE_new_viral_boot.RDS")
estimates_naive <- readRDS("Simulation/SimData/estimates_naive_boot.RDS")
SE_naive <- readRDS("Simulation/SimData/SE_naive_boot.RDS")
estimates_empirical <- readRDS("Simulation/SimData/estimates_empirical_boot.RDS")
SE_empirical <- readRDS("Simulation/SimData/SE_empirical_boot.RDS")

alpha1 = 11
alpha2 = 0.3
alpha3 = 2.8
alpha4 = 1.8
alpha5 = 8
alpha6 = 0.4
alpha7 = 1.3

beta0 <- 8
beta1 <- -2.5



true_value <- c(beta0, beta1)
true_value_viral <- c(alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7)
parnames <- c("beta0","beta1")
parnames_viral <- c("alpha1", "alpha2", "alpha3", "alpha4","alpha5", "alpha6","alpha7")

ncol <- length(parnames)
ncolviral <- length(parnames_viral)

MSE_new <- MSE_empirical <- MSE_naive <- vector("double", ncol)
MSE_ratio_new <- MSE_ratio_empirical <- MSE_ratio_naive <- vector("double", ncol)
SE_sim_new <- SE_sim_empirical <- SE_sim_naive <- vector("double", ncol)
SE_model_new <- SE_model_empirical <- SE_model_naive <- vector("double", ncol)
bias_new <- bias_empirical <- bias_naive <- vector("double", ncol)
bias_ratio_new <- bias_ratio_empirical <- bias_ratio_naive <- vector("double", ncol)
CI_coverage_new <- CI_coverage_empirical <- CI_coverage_naive <- data.frame(matrix(nrow = num_sim, ncol = ncol))
CI_coverage_new_bootstrap_par <- data.frame(matrix(nrow=num_sim, ncol=ncol))
SE_bootstrap_par_avg_new <- vector("double", ncol)

MSE_new_viral <-  MSE_ratio_new_viral <- SE_sim_new_viral <- SE_model_new_viral <- bias_new_viral <- bias_ratio_new_viral <- vector("double", ncolviral)
CI_coverage_new_viral <- data.frame(matrix(nrow = num_sim, ncol = ncolviral))
CI_coverage_new_viral_bootstrap_par <- data.frame(matrix(nrow = num_sim, ncol = ncolviral))
SE_bootstrap_par_avg_new_viral <- vector("double", ncolviral)


CI_lower_new<- estimates_new-1.96*SE_new
CI_upper_new<- estimates_new+1.96*SE_new
CI_lower_naive <- estimates_naive-1.96*SE_naive
CI_upper_naive <- estimates_naive+1.96*SE_naive
CI_lower_empirical <- estimates_empirical-1.96*SE_empirical
CI_upper_empirical <- estimates_empirical+1.96*SE_empirical
CI_lower_new_viral <- estimates_new_viral-1.96*SE_new_viral
CI_upper_new_viral <- estimates_new_viral+1.96*SE_new_viral


CI_lower_new_bootstrap_par <- estimates_new-1.96*bootstrap_par_SE_new
CI_upper_new_bootstrap_par <- estimates_new+1.96*bootstrap_par_SE_new
CI_lower_new_viral_bootstrap_par <- estimates_new_viral-1.96*bootstrap_par_SE_new_viral
CI_upper_new_viral_bootstrap_par <- estimates_new_viral+1.96*bootstrap_par_SE_new_viral

for (i in seq_along(estimates_new)){
  bias_new[i]=sum(estimates_new[,i]-true_value[i])/num_sim
  bias_ratio_new[i]=bias_new[i]/true_value[i]
  SE_sim_new[i] <- sd(estimates_new[,i])
  SE_model_new[i] <- mean(SE_new[,i])
  MSE_new[i]=(bias_new[i])^2+(SE_model_new[i])^2
  MSE_ratio_new[i]=MSE_new[i]/true_value[i]
  SE_bootstrap_par_avg_new[i] <- mean(bootstrap_par_SE_new[,i])
}

for (i in seq_along(estimates_new_viral)){
  bias_new_viral[i]=sum(estimates_new_viral[,i]-true_value_viral[i])/num_sim
  bias_ratio_new_viral[i]=bias_new_viral[i]/true_value_viral[i]
  SE_sim_new_viral[i] <- sd(estimates_new_viral[,i])
  SE_model_new_viral[i] <- mean(SE_new_viral[,i])
  MSE_new_viral[i]=(bias_new_viral[i])^2+(SE_model_new_viral[i])^2
  MSE_ratio_new_viral[i]=MSE_new_viral[i]/true_value_viral[i]
  SE_bootstrap_par_avg_new_viral[i] <- mean(bootstrap_par_SE_new_viral[,i])
}

for (i in 1:num_sim){
  CI_coverage_new[i,] <- (true_value>=CI_lower_new[i,])&(true_value<=CI_upper_new[i,])
  CI_coverage_new_viral[i,] <- (true_value_viral>=CI_lower_new_viral[i,])&(true_value_viral<=CI_upper_new_viral[i,])
  CI_coverage_new_bootstrap_par[i,]<-(true_value>=CI_lower_new_bootstrap_par[i,])&(true_value<=CI_upper_new_bootstrap_par[i,])
  CI_coverage_new_viral_bootstrap_par[i,]<-(true_value_viral>=CI_lower_new_viral_bootstrap_par[i,])&(true_value_viral<=CI_upper_new_viral_bootstrap_par[i,])
}
result_new <- rbind(true_value,sapply(estimates_new,mean),SE_sim_new,SE_model_new,SE_bootstrap_par_avg_new,MSE_new,sqrt(MSE_new),
                    MSE_ratio_new,bias_new,bias_ratio_new,sapply(CI_coverage_new, mean),sapply(CI_coverage_new_bootstrap_par, mean))
result_new_viral <- rbind(true_value_viral,sapply(estimates_new_viral,mean),SE_sim_new_viral,SE_model_new_viral,SE_bootstrap_par_avg_new_viral, MSE_new_viral,sqrt(MSE_new_viral),
                          MSE_ratio_new_viral,bias_new_viral,bias_ratio_new_viral,sapply(CI_coverage_new_viral, mean), sapply(CI_coverage_new_viral_bootstrap_par, mean))

colnames(result_new) <-  parnames
colnames(result_new_viral) <- parnames_viral
rownames(result_new) <- rownames(result_new_viral) <- c("true value","averaged estimates","Simulation SE","model SE","Bootstrap SE","MSE","sqrt MSE","MSE ratio","bias","bias ratio","CI coverage","Bootstrap CI coverage")

for (i in seq_along(estimates_naive)){
  bias_naive[i]=sum(estimates_naive[,i]-true_value[i])/num_sim
  bias_ratio_naive[i]=bias_naive[i]/true_value[i]
  SE_sim_naive[i] <- sd(estimates_naive[,i])
  SE_model_naive[i] <- mean(SE_naive[,i])
  MSE_naive[i]=(bias_naive[i])^2+(SE_model_naive[i])^2
  MSE_ratio_naive[i]=MSE_naive[i]/true_value[i]
}

for (i in 1:num_sim){
  CI_coverage_naive[i,]<-(true_value>=CI_lower_naive[i,])&(true_value<=CI_upper_naive[i,])
}
result_naive <- rbind(true_value,sapply(estimates_naive,mean),SE_sim_naive,SE_model_naive,MSE_naive,sqrt(MSE_naive),MSE_ratio_naive,bias_naive,bias_ratio_naive,sapply(CI_coverage_naive, mean))
rownames(result_naive) = c("true value","averaged estimates","Simulation SE","model SE","MSE","sqrt MSE","MSE ratio","bias","bias ratio","CI coverage")
colnames(result_naive)=parnames

for (i in seq_along(estimates_empirical)){
  bias_empirical[i]=sum(estimates_empirical[,i]-true_value[i])/num_sim
  bias_ratio_empirical[i]=bias_empirical[i]/true_value[i]
  SE_sim_empirical[i] <- sd(estimates_empirical[,i])
  SE_model_empirical[i] <- mean(SE_empirical[,i])
  MSE_empirical[i]=(bias_empirical[i])^2+(SE_model_empirical[i])^2
  MSE_ratio_empirical[i]=MSE_empirical[i]/true_value[i]
}

for (i in 1:num_sim){
  CI_coverage_empirical[i,]<-(true_value>=CI_lower_empirical[i,])&(true_value<=CI_upper_empirical[i,])
}
result_empirical <- rbind(true_value,sapply(estimates_empirical,mean),SE_sim_empirical,SE_model_empirical,MSE_empirical,sqrt(MSE_empirical),MSE_ratio_empirical,bias_empirical,bias_ratio_empirical,sapply(CI_coverage_empirical, mean))
rownames(result_empirical) = c("true value","averaged estimates","Simulation SE","model SE","MSE","sqrt MSE","MSE ratio","bias","bias ratio","CI coverage")
colnames(result_empirical)=parnames 

saveRDS(result_new, "Simulation/SettingI/result_new_B200.RDS")
saveRDS(result_new_viral, "Simulation/SettingI/result_new_viral_B200.RDS")
saveRDS(result_naive, "Simulation/SettingI/result_naive_B200.RDS")
saveRDS(result_empirical, "Simulation/SettingI/result_empirical_B200.RDS")

result_new <- readRDS("Simulation/SettingI/result_new_B200.RDS")
result_new_viral <- readRDS("Simulation/SettingI/result_new_viral_B200.RDS")
result_naive <- readRDS("Simulation/SettingI/result_naive_B200.RDS")
result_empirical <- readRDS("Simulation/SettingI/result_empirical_B200.RDS")

table_names <- c("Method", "Parameter","True Value", "Estimates","SE$_{FIM}$", "SE$_{B}$","SE$_{MC}$",  "Bias (%)","rMSE (%)", "Coverage")
table <- data.frame(matrix(nrow=13, ncol=length(table_names)))
colnames(table) <- table_names
table$Parameter <- c("$\\alpha_1$","$\\alpha_2$","$\\alpha_3$","$\\alpha_4$","$\\alpha_5$","$\\alpha_6$","$\\alpha_7$","$\\beta_0$",NA,NA,"$\\beta_w$",NA,NA)


table$`True Value` <- c(11, 0.3, 2.8, 1.8, 8, 0.4, 1.3, 8,NA,NA,-2.5,NA,NA)
table$Method <- c("NLME",rep(NA,6), rep(c("NLME","Naive","LME"), 2))
table$Estimates <- c(result_new_viral[2,], rbind(result_new[2,], result_naive[2,], result_empirical[2,]))
table$`SE$_{B}$` <- c(result_new_viral[5,], rbind(result_new[5,], rep(NA, 2), rep(NA,2)))
table$`SE$_{FIM}$` <- c(result_new_viral[4,], rbind(result_new[4,], result_naive[4,], result_empirical[4,]))
table$`SE$_{MC}$` <- c(result_new_viral[3,], rbind(result_new[3,], result_naive[3,], result_empirical[3,]))
table$`Bias (%)` <- c(result_new_viral[10,], rbind(result_new[10,], result_naive[9,], result_empirical[9,]))*100
table$`rMSE (%)` <- c(result_new_viral[8,], rbind(result_new[8,], result_naive[7,], result_empirical[7,]))*100
table$`Coverage` <- c(result_new_viral[12,], rbind(result_new[12,], result_naive[10,], result_empirical[10,]))


latex_table<-xtable(table, type = "latex",align=c("ccccccccccc"))
digits(latex_table) <- c(0,0,2,3,3,3,3,3,3,3,2)
print(latex_table, file = "SettingI/latex_table_B200.tex",include.rownames=FALSE,sanitize.text.function = function(x){x},hline.after = c(-1,0,7, 13))



