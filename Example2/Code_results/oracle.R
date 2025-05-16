#install.packages("ALassoSurvIC")
#install.packages("glmnet")
#install.packages("doParallel")
#install.packages("foreach")
#install.packages("Matrix")
#install.packages("devtools")
#install.packages("Metrics")
#install.packages("ggplot2")
#install.packages("cowplot")
#install.packages("tidyr")
#install.packages("robustbase")
#install.packages("readxl")
install.packages("xlsx")
install.packages("BetaBinom")
install.packages("poisson")

rm(list = ls()) 
library(xlsx)
library(readxl)
library(ALassoSurvIC)
library(doParallel)
library(foreach)
library(MASS)
library(glmnet)
library(Matrix)
library(survival)
#library(devtools)
library(Metrics)
library(ggplot2)
library(cowplot)
library(tidyr)
library(robustbase)
library(rmutil)
library(BetaBinom)
library(poisson)

#install_github("feizhe/SSHDI")

setwd("C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/Code/Results/IC30/500x1000_Identity/complete")

#setwd("~/Rui/Gang")
#setwd("F:/R/Rui")

# Loading
beta_n200p200_Ind <- read_excel("Oracle.xlsx", col_names=FALSE)
bet_Oracle <- t(beta_n200p200_Ind)


# Oracle estimation
oracle <- round(colMeans(bet_Oracle), 3)
cat('Oracle_n200p500_Ind', oracle)
cat("\n")



# Empirical SE for oracle
ese_oracle <- c()
for (i in 1:length(bet_Oracle[1, ])) {
  #  ese_oracle[i] <- sd(bet.oracle.all[, i])/sqrt(sum(!is.na(bet.oracle.all[, i])))
  ese_oracle[i] <- sd(bet_Oracle[, i])   
}
cat('Empirical SE for oracle:', round(ese_oracle, 3))
cat("\n")

