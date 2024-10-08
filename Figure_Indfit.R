library(ggpubr)
library(tidyverse)
library(ggplot2)
library(splines)
library(lixoftConnectors)
initializeLixoftConnectors(software="monolix")
setwd(here::here())

### For joint models 
# Import the cleaned data and individual fits
indparest <- read.csv(here::here("DataAnalysis/Results/indest_with_cov_new.csv"), header = TRUE)
indparestemp <- read.csv(here::here("DataAnalysis/Results/indest_with_cov_empirical.csv"), header = TRUE)
data <- read.csv(here::here("Data/Cleaned_data/full2dat_cleaned_with_cov.csv"), header = TRUE)
data_emp <- data %>% 
  mutate(day = day/30)
results_new <- read.csv(here::here("DataAnalysis/Results/results_with_cov_new.csv"), header = TRUE)
results_empirical <- read.csv(here::here("DataAnalysis/Results/results_with_cov_empirical.csv"), header = TRUE)
SE <- results_new$Standard.Error
# patid 610282 
pat1 <- data %>% filter(patid == 610282 & name == "rna") 
indparest1 <- indparest %>% filter(id == 610282)
P1 <- indparest1$p1
lambda1 <- indparest1$b1
gamma <- results_new$Estimates[3]
P2 <- indparest1$p2
beta1 <- indparest1$beta1
beta2 <- indparest1$beta2
beta3 <- indparest1$beta3
beta4 <- indparest1$beta4
age <- unique(pat1$age_st)
decaylast <- unique(pat1$Ti)
rebound1st <- min(pat1[which(pat1$decay == 0), ]$day)
datamid <- data.frame(day = seq(decaylast, rebound1st, 0.1))
datamid <- datamid %>% 
  mutate(ydecay = log10( exp(P1-(lambda1+gamma * age)*day) + exp(P2))) %>% 
  mutate(yrebound = beta1 * (day-decaylast) / ((day-decaylast)  + exp(beta2-beta3*(day-decaylast))) + beta4) %>% 
  group_by(day) %>% 
  mutate(observation = max(ydecay, yrebound))

fun1 <- function(x){
  ifelse(x<=decaylast, log10( exp(P1-(lambda1+gamma*age)*x) + exp(P2)), NA)
}
fun2 <- function(x){
  ifelse(x>=rebound1st, beta1 * (x-decaylast) / ((x-decaylast)  + exp(beta2-beta3*(x-decaylast))) + beta4, NA)
}
indparestemp1 <- indparestemp %>% filter(id == 610282)
alpha0 <- indparestemp1$alpha0
alpha1 <- indparestemp1$alpha1
alpha2 <- indparestemp1$alpha2
gamma.emp <- results_empirical$Estimates[2]
fun3 <- function(x){
  alpha0 + gamma.emp * age + alpha1 * x/30 + alpha2 * (x/30)^2
}

(p1 <- ggplot(pat1, aes(x = day, y = observation)) + geom_point(size=3) + 
    stat_function(fun = fun1, aes(linetype = "Joint")) + stat_function(fun = fun2) + geom_smooth(data=datamid, se=FALSE, color = "black",size=0.5) +
    stat_function(fun = fun3, aes(linetype = "Empirical Model")) +
    scale_linetype_manual(values = c("dashed","solid"), name = "Methods", guide = FALSE) +
    annotate("text", x = 150, y=1.8, label = "detection limit", size =5) +
    xlab("Days") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"), limits = c(0,6)) + theme_classic() + theme(text = element_text(size = 20)) + 
    geom_hline(yintercept=2, linetype="dotted")
)  


# patid 271775
ind <- 271775
pat3 <- data %>% filter(patid == ind & name == "rna") %>% 
  mutate(observation = ifelse(censor == 1, log10(50), observation))
indparest3 <- indparest %>% filter(id == ind)
P1 <- indparest3$p1
lambda1 <- indparest3$b1
P2 <- indparest3$p2
beta1 <- indparest3$beta1
beta2 <- indparest3$beta2
beta3 <- indparest3$beta3
beta4 <- indparest3$beta4
age <- unique(pat3$age_st)
decaylast <- unique(pat3$Ti)
ind3t56 <- log10(exp(P1-(lambda1+gamma*age)*56) + exp(P2))
ind3t84 <-  beta1 * (84-56) / (84-56  + exp(beta2-beta3*(84-56))) + beta4
rebound1st <- min(pat3[which(pat3$decay == 0), ]$day)
datamid <- data.frame(day = seq(decaylast-1, rebound1st+1, 0.1))
datamid <- datamid %>% 
  mutate(ydecay = log10( exp(P1-(lambda1+gamma*age)*day) + exp(P2))) %>% 
  mutate(yrebound = beta1 * (day-decaylast) / ((day-decaylast)  + exp(beta2-beta3*(day-decaylast))) + beta4) %>% 
  mutate(observation = ifelse(day <= 0.5*(decaylast+rebound1st), ydecay, yrebound))

fun1 <- function(x){
  ifelse(x<=decaylast, log10( exp(P1-(lambda1+gamma*age)*x) + exp(P2)), NA)
}
fun2 <- function(x){
  ifelse(x>=rebound1st, beta1 * (x-decaylast) / ((x-decaylast)  + exp(beta2-beta3*(x-decaylast))) + beta4, NA)
}


indparestemp3 <- indparestemp %>% filter(id == ind)
alpha0 <- indparestemp3$alpha0
alpha1 <- indparestemp3$alpha1
alpha2 <- indparestemp3$alpha2

fun3 <- function(x){
  alpha0 + gamma.emp * age + alpha1 * x/30 + alpha2 * (x/30)^2
}

(p3 <- ggplot(pat3, aes(x = day, y = observation)) + 
    geom_point(aes(shape=factor(censor)),size=3) + 
    scale_shape_manual(name = "Censored", values = c(16, 1), labels=c('No', 'Yes'))+
    stat_function(fun = fun1, aes(linetype = "NLME model")) + stat_function(fun = fun2) + geom_smooth(data=datamid, se=FALSE, color = "black",size=0.5, span = 1) +
    stat_function(fun = fun3, aes(linetype = "LME model")) +
    annotate("text", x = 158, y=1.8, label = "detection limit", size = 5) +
    scale_linetype_manual(values = c("dashed","solid"), name = "Methods") +
    geom_hline(yintercept=2, linetype="dotted") + 
    xlab("Days") + scale_y_continuous(bquote("Viral load (in" ~ log[10]~"-scale)"), limits = c(0, 5)) + theme_classic()  + theme(text = element_text(size = 20)))

# Nonparametric bootstrap to find prediction interval
allpatid <- unique(data$patid)
n <- length(allpatid)
pat3.ind <- data.frame()
for (B in 1:200){
  print(glue::glue("B = ",B))
  datanew <- data.frame()
  for (i in 1:n){
    sampled_patient <- sample(allpatid, 1)
    # Create a new data frame with the sampled patient data
    sampled_data <- data[which(data$patid == sampled_patient), ] %>% 
      mutate(ID = paste0(sampled_patient,".",i))
    datanew <- rbind(datanew, sampled_data)
  }
  write.csv(datanew, "Data/Other/sampled_data.csv", row.names = FALSE)
  # Read in data for Monolix
  data_new = list(dataFile = paste0('Data/Other/sampled_data.csv'),
                  headerTypes =c("ignore","time","ignore","cens","limit","regressor","regressor","ignore", "ignore","ignore", "contcov","contcov","contcov","contcov", "obsid","observation","regressor","id"),
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
  setCorrelationBlocks(ID = list(c("alpha2", "alpha4", "alpha5", "alpha6")))
  setIndividualParameterVariability(alpha1 = TRUE, alpha2 = TRUE, alpha3 = TRUE, alpha4 = TRUE, alpha5 = FALSE, alpha6 = TRUE, alpha7 = FALSE, beta0 = TRUE, beta1 = FALSE)
  setIndividualParameterDistribution(alpha1="normal",alpha2="normal",alpha3="normal",alpha4="normal",alpha5="normal",alpha6="normal",alpha7="normal", beta0="normal", beta1="normal")
  setPopulationParameterInformation(alpha1_pop = list(initialValue = 11), 
                                    alpha2_pop = list(initialValue = 0.3),
                                    alpha3_pop = list(initialValue = 5.7),
                                    alpha4_pop = list(initialValue = 1.4),
                                    alpha5_pop = list(initialValue = 5),
                                    alpha6_pop = list(initialValue = 0.2),
                                    alpha7_pop = list(initialValue = 2.4),
                                    beta0_pop = list(initialValue = 9.7), 
                                    beta1_pop  = list(initialValue = -2))
  setCovariateModel(alpha2 = c(age_st = TRUE))
  # run the estimation
  setScenario(scenario)
  runScenario()
  # get individual parameters for plots later
  indest.boot <- getEstimatedIndividualParameters()$conditionalMean
  indest.boot <- indest.boot[grep("271775", indest.boot$id), ] %>% 
    mutate(gamma = getEstimatedPopulationParameters()[3])
  pat3.ind <- rbind(pat3.ind, indest.boot)
}
saveRDS(pat3.ind, "Data/Other/pat3ind.RDS")
pat3.ind <- readRDS("Data/Other/pat3ind.RDS")

pat3.ind <- pat3.ind %>% 
  mutate(t56 = log10(exp(p1-(b1+gamma*age)*56) + exp(p2))) %>% 
  mutate(t84 = beta1 * (84-56) / (84-56  + exp(beta2-beta3*(84-56))) + beta4)
lower.56 <- ind3t56 - 1.96 * sd(pat3.ind$t56)
upper.56 <- ind3t56 + 1.96 * sd(pat3.ind$t56)
lower.84 <- ind3t84 - 1.96 * sd(pat3.ind$t84)
upper.84 <- ind3t84 + 1.96 * sd(pat3.ind$t84)

p3 + geom_segment(aes(x = 56, y = lower.56, xend = 56, yend = upper.56))+
  geom_segment(aes(x = 54, y = lower.56, xend = 58, yend = lower.56))+
  geom_segment(aes(x = 54, y = upper.56, xend = 58, yend = upper.56))+
  geom_point(aes(x=56, y=ind3t56), shape = 4, size = 4) + 
  geom_point(aes(x=84, y=ind3t84), shape = 4, size = 4)+ 
  geom_segment(aes(x = 84, y = lower.84, xend = 84, yend = upper.84))+
  geom_segment(aes(x = 82, y = lower.84, xend = 86, yend = lower.84))+
  geom_segment(aes(x = 82, y = upper.84, xend = 86, yend = upper.84))
