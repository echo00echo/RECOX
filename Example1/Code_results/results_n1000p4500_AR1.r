install.packages("openxlsx")
install.packages("magrittr")

#library(readxl)
library(openxlsx)
library(magrittr)
library("readxl")
library(plotrix)
#library(readxl)
library(openxlsx)
library(magrittr)
#library(readxl)
library(openxlsx)
library(magrittr)
library(readxl)
library(ggplot2)
library(openxlsx)
library(magrittr)
library(tikzDevice)
library(latex2exp)
library(purrr)
library(tibble)
library(cowplot)
library(installr)
library(qqplotr)


rm(list=ls())
setwd("C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/Code/Results/Case3_simulation/n1000p4500/Hera/20230306_1_46_100")


#a <- loadWorkbook('C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/Code/1026/yourExcel1.xlsx')
#b <- loadWorkbook('C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/Code/1026/yourExcel5.xlsx')
a <- loadWorkbook('C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/Code/Results/Case3_simulation/n1000p4500/Hera/20230306_1_46_100/n1000p4500_beta_p_ar1.xlsx')
b <- loadWorkbook('C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/Code/Results/Case3_simulation/n1000p4500/Hera/20230306_1_46_100/n1000p4500_sample_p_ar1.xlsx')
c <- loadWorkbook('C:/Users/echo0/Box Sync/RuiYang/Project3-IntervalCensoring/code/Results/Case3_simulation/n1000p4500/Hera/20230306_1_46_100/n1000p4500_screening_sample_p.xlsx')
#sheetNames <- sheets(a)
#for(i in 1:length(sheetNames))
#{
#  assign(sheetNames[i],readWorkbook(a,sheet = i, colNames = FALSE, skipEmptyRows = FALSE,
#                                    skipEmptyCols = FALSE))
#}

MSE.LASSO <- TP.LASSO_u <- TP.LASSO_b <- FP.LASSO_u <- FP.LASSO_b <- c()
CP_u <- CP_b <- alpha_u <-alpha_b <- LL_u <- LL_b <- UL_u <- UL_b <- Bias <- sel.freq.all <- bet.oracle.all <-  set <- set.rec <- bet.rec <- c()
pvalues_u <- pvalues_b <- qvalues_u <- qvalues_b <- ase_u <- ase_b <- ese <- ese_oracle <- vars_u.all <- vars_b.all <- c()
sel.freq.all_b_qv <- sel.freq.all_u_qv <- sel.freq.all_b_pv <- sel.freq.all_u_pv <- c()
nzero_u_pv.all <- nzero_u_qv.all <- nzero_b_pv.all<- nzero_b_qv.all <- c()
qvalues_final_u <- qvalues_final_b <- c() 


# beta = [-1.2; zeros(48,1); -0.8; zeros(49,1); 0.6; zeros(49,1); 0.9; zeros(49,1); 1.5];

for (k in 1:100){
  bet <- bet0 <- c(-1.2, rep(0,48), -0.8, rep(0,49), 0.6, rep(0,49), 0.9, rep(0,49), 1.5, rep(0, 800)) # True value beta0
  # c(rep(0.5,2),rep(0,6),rep(0.5,2)) 
  
  s0 <- which(bet!=0) # True active set
  p <- length(bet) #dim of Beta
  Beta_Hat.all1 <- readWorkbook(a, sheet = k, colNames = FALSE, skipEmptyRows = FALSE,
                                skipEmptyCols = FALSE)
  

  Ycount <- readWorkbook(b, sheet = k, colNames = FALSE, skipEmptyRows = FALSE,
                         skipEmptyCols = FALSE)

  Bet_SPARES <- rowMeans(Beta_Hat.all1)
  
  coef2 <- readWorkbook(c, sheet = k, colNames = FALSE, skipEmptyRows = FALSE,
                      skipEmptyCols = FALSE)
  
  set.rec <- c()
  
  for (t in 1:400) {
    set <- which(coef2[, t]!=0)
    set.rec <- c(set.rec, unlist(set))
  }
 
#  set<-(1:p)[coef2!=0]
#  set.rec <- c(set.rec, unlist(set)) #Store set from each SPARE step 
  
  
  # estimator coefficient variance
  est.var_u <- function(beta,Ycount, unbiased =TRUE){
    n = dim(Ycount)[1]
    B = dim(Ycount)[2]
    n1 = sum(Ycount[,1])
    
    var1 <- (n/(n-n1))^2 * sum(cov(beta,t(Ycount))^2)
    var2 <- ifelse(unbiased, var1 - n1/(n-n1)*n/B*var(beta), var1)
    return(var2)
  }
  
  
  est.var_b <- function(beta,Ycount, unbiased =FALSE){
    n = dim(Ycount)[1]
    B = dim(Ycount)[2]
    n1 = sum(Ycount[,1])
    
    var1 <- (n/(n-n1))^2 * sum(cov(beta,t(Ycount))^2)
    #  var2 <- ifelse(unbiased, var1 - n1/(n-n1)*n/B*var(beta), var1)
    return(var1)
  }
  
  
  ### The p-value of testing H0: B_j*=0
  vars_b <-apply(Beta_Hat.all1, 1, est.var_b, Ycount=(Ycount)) # biased variance estimator for j=1,...,p
  
  vars_u <-apply(Beta_Hat.all1, 1, est.var_u, Ycount=(Ycount)) # biased variance estimator for j=1,...,p
  
  
  # biased Var from each replication
  vars_b.all <-  rbind(vars_b.all, c(vars_b)) 
  
  
  # unbiased Var from each replication
  vars_u.all <-  rbind(vars_u.all, c(vars_u)) 
  
  
  
  sds_u <- sqrt(vars_u)
  cat("Check negative SD for each replication:", sds_u)
  
  
  sds_b <- sqrt(vars_b)
  
# cat("\n")
#  if(sum(is.na(sds_u))>0) next # if negative SD then skip to the next k
  
  
  
  pvs_u<-2*(1-pnorm(abs(Bet_SPARES)/sds_u))  #p-value unbiased
  nzero_u_pv <- (1:p)[pvs_u<0.05/p & !is.na(pvs_u)] #non-zero signals chosen by unbiased p-values and Bonferroni correction

  nzero_u_pv.all <-  c(nzero_u_pv.all, unlist(nzero_u_pv))#non-zero signals chosen by unbiased p-values and Bonferroni correction at each replication
  
  
  
  
  
  # Selection frequency (from testing) for each beta_j based on unbiased pvalue 
  temptable <- table(nzero_u_pv.all) # frequency of each beta.j to be selected
  sel_freq_u_pv<-rep(0,p)
  #B = dim(Ycount)[2]
  nloop = k
  sel_freq_u_pv[as.numeric(names(temptable))]<-temptable/nloop
 # sel.freq.all_u_pv <-  rbind(sel.freq.all_u_pv, c(sel_freq_u_pv)) # selection frequency at each MC step
  s_freq_u_pv <- round(sel_freq_u_pv,3)
  cat('Selection frequency from testing_unbiased_pv:', s_freq_u_pv[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_s_freq_u_pv <- round(mean(s_freq_u_pv[c(2:49, 51:99, 101:149, 151:199, 201:1000)]), 3) 
  cat('# Selection frequency for 0s from testing_unbiased_pv:', ave_0_s_freq_u_pv)
  cat("\n")
  
  
  
  
  
  
  
  
  
  pvs_b<-2*(1-pnorm(abs(Bet_SPARES)/sds_b))  #p-value biased
  nzero_b_pv <- (1:p)[pvs_b<0.05/p & !is.na(pvs_b)] #non-zero signals chosen by biased p-values and Bonferroni correction

  
  
  nzero_b_pv.all <-  c(nzero_b_pv.all, unlist(nzero_b_pv))#non-zero signals chosen by biased p-values and Bonferroni correction at each replication
 
  
  
  # Selection frequency (from testing) for each beta_j based on biased pvalue 
  temptable <- table(nzero_b_pv.all) # frequency of each beta.j to be selected
  sel_freq_b_pv<-rep(0,p)
#  B = dim(Ycount)[2]
  nloop = k
  sel_freq_b_pv[as.numeric(names(temptable))]<-temptable/nloop
 # sel.freq.all_b_pv <-  rbind(sel.freq.all_b_pv, c(sel_freq_b_pv)) # selection frequency at each MC step
  s_freq_b_pv <- round(sel_freq_b_pv,3)
  cat('Selection frequency from testing_biased_pv:', s_freq_b_pv[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_s_freq_b_pv <- round(mean(s_freq_b_pv[c(2:49, 51:99, 101:149, 151:199, 201:1000)]), 3) 
  cat('# Selection frequency for 0s from testing_biased_pv:', ave_0_s_freq_b_pv)
  cat("\n")
  
  
  
  


  
  
  
  # unbiased pvalues from each replication
  pvalues_u <-  rbind(pvalues_u, c(pvs_u)) 
  
  
  # biased pvalues from each replication
  pvalues_b <-  rbind(pvalues_b, c(pvs_b)) 
  

  
  # active set from each replication - unbiased - pvalue
  nzero_u_pv #non-zero signals chosen by unbiased p-values and Bonferroni correction
  cat("Active set from each replication-unbiased-pvalue:", nzero_u_pv)
  
  cat("\n")
  
  
  # active set from each replication - biased -pvalue
  nzero_b_pv  #non-zero signals chosen by biased p-values and Bonferroni correction
  cat("Active set from each replication-biased:", nzero_b_pv)
  
  cat("\n")
  

  
  
  ### unbiased TP & FP based on the test
  #TP.LASSO_u <- c(TP.LASSO_u, sum(as.numeric(nzero_u %in% c(1, 2, 6, 7, 10))))
  #FP.LASSO_u <- c(FP.LASSO_u, sum(as.numeric(nzero_u %in% c(3:5, 8, 9))))
  #cat('TP_u:', mean(TP.LASSO_u))  
  #cat("\n")
  #cat('FP_u:', mean(FP.LASSO_u))  
  #cat("\n")
  
  
  ### biased TP & FP based on the test
  #TP.LASSO_b <- c(TP.LASSO_b, sum(as.numeric(nzero_b %in% c(1, 2, 6, 7, 10))))
  #FP.LASSO_b <- c(FP.LASSO_b, sum(as.numeric(nzero_b %in% c(3:5, 8, 9))))
  #cat('TP_b:', mean(TP.LASSO_b))  
  #cat("\n")
  #cat('FP_b:', mean(FP.LASSO_b))  
  #cat("\n")
  
  
  
  # asympotic SE from each replication
  ase_u <-  rbind(ase_u, c(sds_u)) 
  
  
  
  # asympotic SE from each replication
  ase_b <-  rbind(ase_b, c(sds_b)) 
  
  
  # mean of asympotic SE - unbiased
  ase_final_u <- round(colMeans(ase_u, na.rm=TRUE), 3)
  cat('Average of asympotic SE - unbiased:', ase_final_u[c(1, 50, 100, 150, 200)])
  cat("\n")
  
  
  ave_0_ase_final_u <- round(mean(ase_final_u[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3) 
  cat('Average of asympotic SE for 0s - unbiased:', ave_0_ase_final_u)
  cat("\n")
  
  
  # mean of asympotic SE - biased
  ase_final_b <- round(colMeans(ase_b, na.rm=TRUE), 3)
  cat('Average of asympotic SE - biased:', ase_final_b[c(1, 50, 100, 150, 200)])
  cat("\n")
  
  
  ave_0_ase_final_b <- round(mean(ase_final_b[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3) 
  cat('Average of asympotic SE for 0s - biased:', ave_0_ase_final_b)
  cat("\n")
  
  
  
  
  
  # asymptotic 95% CI for each beta.j* - unbiased
  z <- qnorm(0.975,mean=0,sd=1)
  
  ll_u <- Bet_SPARES- z%*%sds_u
  ul_u <- Bet_SPARES + z%*%sds_u
  
  LL_u <- rbind(LL_u, c(ll_u))
  UL_u <- rbind(UL_u, c(ul_u))
  
  L_u <- round(colMeans(LL_u)[c(1,50,100,150,200)], 3)
  U_u <- round(colMeans(UL_u)[c(1,50,100,150,200)], 3)
  ave_0_L_u <- mean(colMeans(LL_u, na.rm=TRUE)[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE)
  ave_0_U_u <- mean(colMeans(UL_u, na.rm=TRUE)[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE) 
  
  
  
  # asymptotic 95% CI for each beta.j* - biased
  z <- qnorm(0.975,mean=0,sd=1)
  
  ll_b <- Bet_SPARES- z%*%sds_b
  ul_b <- Bet_SPARES + z%*%sds_b
  
  LL_b <- rbind(LL_b, c(ll_b))
  UL_b <- rbind(UL_b, c(ul_b))
  
  
  
  # Selection frequency
  temptab<-table(set.rec) # frequency of each beta.j to be selected from the screening step
  sel_freq<-rep(0,p)
  B = dim(Ycount)[2]
  sel_freq[as.numeric(names(temptab))]<-temptab/B
  sel.freq.all <-  rbind(sel.freq.all, c(sel_freq)) # selection frequency at each MC step
  s_freq <- round(colMeans(sel.freq.all),3)
  cat('Selection frequency:', s_freq[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_s_freq <- round(mean(s_freq[c(2:49, 51:99, 101:149, 151:199, 201:1000)]), 3) 
  cat('# Selection frequency for 0s:', ave_0_s_freq)
  cat("\n")
  
  
  # Store Beta hat from each MC step
  bet.rec <- rbind(bet.rec, c(Bet_SPARES)) 
  
  
  
  # Empirical SE for beta hat
  for (i in 1:length(bet.rec[1, ])) {
    #  ese[i] <- sd(bet.rec[, i])/sqrt(sum(!is.na(bet.rec[, i])))
    ese[i] <- sd(bet.rec[, i]) 
  }
  cat('Empirical SE for beta hat:', round(ese,3)[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_ese <- round(mean(ese[c(2:49, 51:99, 101:149, 151:199, 201:1000)]), 3) 
  cat('Empirical SE for beta hat for 0s:', ave_0_ese)
  cat("\n")
  
  
  
  # Final estimator beta hat
  bet.final <- round(colMeans(bet.rec), 3)
  #bet.final[abs(bet.final) < 1e-1] <- 0 #Set beta=0 if <0.0001 (removed! estimations do not need to be sparsed)
  cat('Final estimator beta hat:', bet.final[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_bet.final <- round(mean(bet.final[c(2:49, 51:99, 101:149, 151:199, 201:1000)]), 3) 
  cat('Final estimator beta hat for 0s:', ave_0_bet.final)
  cat("\n")
  
  
  # Coverage Prob (should be 95%) - unbiased
  for (i in 1:length(LL_u[1, ])) {
    CP_u[i] <- sum(LL_u[, i] < bet.final[i] & UL_u[, i]> bet.final[i], na.rm=TRUE)/length(LL_u[, 1])
    
  }
  cat('Coverage Prob - unbiased:', round(CP_u,3)[c(1,50,100,150,200)])
  cat("\n")
  
 # CP_u[i] <- sum(LL_u[, 2] < bet.final[2] & UL_u[, 2]> bet.final[2], na.rm=TRUE)/length(LL_u[, 1])
  
  
  ave_0_cp_u <- round(mean(CP_u[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3) 
  cat('Coverage Prob for 0s - unbiased:', ave_0_cp_u)
  cat("\n")
  
  
  
  # Coverage Prob (should be 95%) - biased
  for (i in 1:length(LL_b[1, ])) {
    CP_b[i] <- sum(LL_b[, i] < bet.final[i] & UL_b[, i]> bet.final[i])/length(LL_b[, 1])
    
  }
  cat('Coverage Prob - biased:', round(CP_b,3)[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_cp_b <- round(mean(CP_b[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3) 
  cat('Coverage Prob for 0s - biased:', ave_0_cp_b)
  cat("\n")
  
  
  
  
  # alpha level (should be 5%) - unbiased: prob of rejecting the null hypothesis when it was true. 
  # i.e. 5 times out of 100, the CI will not contain the estimate
  for (i in 1:length(LL_u[1, ])) {
    alpha_u[i] <- sum(LL_u[, i] > bet.final[i] | UL_u[, i]< bet.final[i], na.rm=TRUE)/length(LL_u[, 1])
    
  }
  
  cat('Alpha Level - unbiased:', round(alpha_u,3)[c(1,50,100,150,200)])
  cat("\n")
  
  # CP_u[i] <- sum(LL_u[, 2] < bet.final[2] & UL_u[, 2]> bet.final[2], na.rm=TRUE)/length(LL_u[, 1])
  
  
  ave_0_alpha_u <- round(mean(alpha_u[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3) 
  cat('alpha level for 0s - unbiased:', ave_0_alpha_u)
  cat("\n")
  
  
  
  # alpha level (should be 5%) - biased: prob of rejecting the null hypothesis when it was true. 
  # i.e. 5 times out of 100, the CI will not contain the estimate
  for (i in 1:length(LL_b[1, ])) {
    alpha_b[i] <- sum(LL_b[, i] > bet.final[i] | UL_b[, i]< bet.final[i])/length(LL_b[, 1])
    
  }
  cat('Alpha Level - biased:', round(alpha_b,3)[c(1,50,100,150,200)])
  cat("\n")
  
  
  ave_0_alpha_b <- round(mean(alpha_b[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3) 
  cat('alpha level for 0s - biased:', ave_0_alpha_b)
  cat("\n")
  
  
  
  # Bias 
  bias <- round(Bet_SPARES - bet0, 3)
  Bias <- rbind(Bias, c(bias)) # bias at each step
  Bias_final <- round(colMeans(Bias),3)
  cat('Bias', Bias_final[c(1,50,100,150,200)])
  
  cat("\n")
  
  
  ave_0_bias <- round(mean(Bias_final[c(2:49, 51:99, 101:149, 151:199, 201:1000)]), 3) 
  cat('Bias for 0s:', ave_0_bias)
  cat("\n")
  
  
  
  # Empirical Type I error rate for 0s
  ave_0_alpha_level <- sum(pvalues_u[, c(2:49, 51:99, 101:149, 151:199, 201:1000)]<=0.05, na.rm=TRUE)/(995*100)
  cat('Empirical Type I error rate for 0s:', ave_0_alpha_level)
  cat("\n")
  
  # Oracle estimation
#  oracle <- round(colMeans(bet.oracle.all), 3)
#  cat('Oracle', oracle[c(1,2,6,7,10)])
#  cat("\n")
  
  
  
#  ave_0_oracle <- round(mean(oracle[c(3:5, 8 ,9)]), 3) 
#  cat('Oracle for 0s:', ave_0_oracle)
#  cat("\n")
  
  
  
  # Empirical SE for oracle
#  for (i in 1:length(bet.oracle.all[1, ])) {
    #  ese_oracle[i] <- sd(bet.oracle.all[, i])/sqrt(sum(!is.na(bet.oracle.all[, i])))
#    ese_oracle[i] <- sd(bet.oracle.all[, i])   
#  }
#  cat('Empirical SE for oracle:', round(ese_oracle, 3)[c(1,2,6,7,10)])
#  cat("\n")
  
  
  
#  ave_0_ese_oracle <- round(mean(ese_oracle[c(3:5, 8 ,9)]), 3) 
#  cat('Empirical SE for oracle:', ave_0_ese_oracle)
#  cat("\n")
  
  
  
  # MSE
#  MSE.LASSO <- c(MSE.LASSO, matrix(c(Bet_SPARES-bet0),1,p)%*%VarCovS%*%t(matrix(c(Bet_SPARES-bet0),1,p)))
#  cat('MSE', round(median(MSE.LASSO),3)) 
#  cat("\n")
#  cat('SD of MSE', round(sd(MSE.LASSO),3)) 
#  cat("\n")
  
  print(k)  
  flush.console()
  
#  if (length(bet.rec[, 1])==100) break
  
} 



# Outputs
#bet.oracle.all
#bet.rec
#pvalues
#ase
#sds
#vars.all

write.csv(vars_b.all,'Vars_b_ar1_B400_n1000.csv')
write.csv(vars_u.all,'Vars_u_ar1_B400_n1000.csv')
write.csv(pvalues_u,'pvalues_u_ar1_B400_n1000.csv')
#write.csv(qvalues_final_u,'qvalues_u_ar1_n1000.csv')
write.csv(pvalues_b,'pvalues_b_ar1_n1000.csv')
#write.csv(qvalues_final_b,'qvalues_b_ar1_n1000.csv')
write.csv(bet.rec,'beta_AR(1)_n1000.csv')
write.csv(ase_u, 'asymptotic_SE_u_ar1_n1000.csv')
write.csv(ase_b, 'asymptotic_SE_b_ar1_n1000.csv')
#write.csv(bet.oracle.all, 'Oracle_ar1_n1000.csv')

mylist <- list("bet"=bet.final, "bet_0s"=ave_0_bet.final,
               "Bias"=round(colMeans(Bias),3), "Bias_0s"=ave_0_bias, 
               "E_SE"=round(ese,3), "E_SE_0s"=ave_0_ese,
               "mean_a_SE_u"=ase_final_u, "mean_a_SE_0s_u"=ave_0_ase_final_u, 
               "mean_a_SE_b"=ase_final_b, "mean_a_SE_0s_b"=ave_0_ase_final_b, 
               "Coverage_Prob_u"=round(CP_u,3), "Coverage_Prob_0s_u"=ave_0_cp_u, 
               "Coverage_Prob_b"=round(CP_b,3), "Coverage_Prob_0s_b"=ave_0_cp_b,
               "sel_freq_u_pv"=s_freq_u_pv[c(1,50, 100, 150,200)], "sel_freq_0s_u_pv"=ave_0_s_freq_u_pv, 
               "sel_freq_b_pv"=s_freq_b_pv[c(1,50, 100, 150, 200)], "sel_freq_0s_b_pv"=ave_0_s_freq_b_pv, 
               "sel_freq"=s_freq[c(1,50,100,150,200)], "sel_freq_0s"=ave_0_s_freq,
               "Alpha Level-unbiased"=round(alpha_u,3)[c(1,50, 100, 150, 200)],
               "Alpha Level-unbiased_0s"=round(mean(alpha_u[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3),
               "Alpha Level-biased"=round(alpha_b,3)[c(1,50, 100, 150, 200)],
               "Alpha Level-biased_0s"=round(mean(alpha_b[c(2:49, 51:99, 101:149, 151:199, 201:1000)], na.rm=TRUE), 3),
               "Empirical Type I error rate"=round(sum(pvalues_u[, c(2:49, 51:99, 101:149, 151:199, 201:1000)]<=0.05 ,na.rm=TRUE)/(995*100), 3))


cat(capture.output(print(mylist), file="ar1_n1000_p4500_IC1_AR1.txt"))

#order(colMeans(pvalues_b))
#round(colMeans(pvalues_b),7)

# 
# 
# est <- bet.final[c(1, 50, 100, 150, 200)]
# truth <- bet[c(1, 50, 100, 150, 200)]
# 
# L_u <- round(colMeans(LL_u)[c(1,50,100,150,200)], 3)
# U_u <- round(colMeans(UL_u)[c(1,50,100,150,200)], 3)
# ave_0_L_u <- mean(colMeans(LL_u, na.rm=TRUE)[c(2:49, 51:99, 101:149, 151:199, 201:4500)], na.rm=TRUE)
# ave_0_U_u <- mean(colMeans(UL_u, na.rm=TRUE)[c(2:49, 51:99, 101:149, 151:199, 201:4500)], na.rm=TRUE) 
# 
# 
# 
# 
# testplot <- read_excel("plot_data.xlsx")
# testplot1  <- read_excel("plot_data_AR1.xlsx")
# testplot2  <- read_excel("plot_data_CS.xlsx")
# 
# 
# plot(x = testplot2$position, y = testplot2$truth, ylim = c(-1.5, 1.8), type = 'p', col = 'red', pch = 1, lwd = 2, cex = 2,
#      xlab = "Position", ylab = "Value", xaxt = "n", main = "n=500, p=1000, CS")
# er <- (testplot2$CIU - testplot2$CIL)/2
# plotCI(x = testplot2$position, y = testplot2$truth, er, err="y", slty=par("lty"), add=TRUE)
# points(x= testplot2$position, y =  testplot2$b_est, type = 'p', col = 'green', pch = 3, lwd = 2, cex = 2)
# abline(h=0, col="blue", lty=2)
# axis(1,                         # Define x-axis manually
#      at = c(0, 50,100, 150, 200, 250),
#      labels = c("1", "50", "100", "150", "200", "averaged 0s'"))
# 
# legend("bottomright", pch = c(1, 3), lty = c(NA, NA), legend = c("Truth", "Estimate"),
#        lwd = 3, col = c("red", "green", "blue"))
# 
# 
# 
# # Make the window wider than taller
# #windows(width = 4.5, height = 4)
# 
# 
# par()
# windows()
# par(mfrow=c(3,1))
# 
# # Save current graphical parameters
# #opar <- par(no.readonly = TRUE)
# 
# 
# 
# # Change the margins of the plot (the fourth is the right margin)
# #par(mar = c(6, 4.1, 4.1, 2.1))
# 
# 
# plot(x = testplot$position, y = testplot$truth, ylim = c(-1.5, 1.8), type = 'b', col = 'red', pch = 16, lwd = 6, xlab = " ", ylab = "Value", xaxt = "n", main = "n=500, p=1000, Identity", cex.main = 3)
# #polygon(c(testplot$position,rev(testplot$position)),c(testplot$CIL,rev(testplot$CIU)),col = rgb(1,1,1), border = NA)
# #plot(x = x_raw, y = y,type = 'p')
# lines(x= testplot$position, y = testplot$CIL, lty = 'dashed', col = 'blue', lwd = 4)
# lines(x= testplot$position, y =  testplot$b_est, type = 'b', col = 'green', pch = 1, lwd = 2)
# lines(x= testplot$position, y =  testplot$CIU, lty = 'dashed', col = 'blue', lwd = 4)
# axis(1,                         # Define x-axis manually
#      at = c(0, 50,100, 150, 200, 250),
#      labels = c("1", "50", "100", "150", "200", "averaged 0s'"))
# # Adding a legend
# #legend("bottomright", lty = c(1, 1, 2), legend = c("Truth", "Estimate", "95% CI"),
# #       lwd = 3, col = c("red", "green", "blue"))
# 
# 
# 
# 
# 
# 
# plot(x = testplot1$position, y = testplot1$truth, ylim = c(-1.5, 1.8), type = 'b', col = 'red', pch = 16, lwd = 6, xlab = " ", ylab = "Value", xaxt = "n", main = "n=500, p=1000, AR(1)", cex.main = 3)
# #polygon(c(testplot1$position,rev(testplot1$position)),c(testplot1$CIL,rev(testplot1$CIU)),col = rgb(1,1,1), border = NA)
# #plot(x = x_raw, y = y,type = 'p')
# lines(x= testplot1$position, y = testplot1$CIL, lty = 'dashed', col = 'blue', lwd = 4)
# lines(x= testplot1$position, y =  testplot1$b_est, type = 'b', col = 'green', pch = 1, lwd = 2)
# lines(x= testplot1$position, y =  testplot1$CIU, lty = 'dashed', col = 'blue', lwd = 4)
# axis(1,                         # Define x-axis manually
#      at = c(0, 50,100, 150, 200, 250),
#      labels = c("1", "50", "100", "150", "200", "averaged 0s'"))
# # Adding a legend
# #legend("bottomright", lty = c(1, 1, 2), legend = c("Truth", "Estimate", "95% CI"),
# #       lwd = 3, col = c("red", "green", "blue"))
# 
# 
# 
# plot(x = testplot2$position, y = testplot2$truth, ylim = c(-1.5, 1.8), type = 'b', col = 'red', pch = 16, lwd = 6, xlab = " ", ylab = "Value", xaxt = "n", main = "n=500, p=1000, CS", cex.main = 3)
# #polygon(c(testplot2$position,rev(testplot2$position)),c(testplot2$CIL,rev(testplot2$CIU)),col = rgb(1,1,1), border = NA)
# #plot(x = x_raw, y = y,type = 'p')
# lines(x= testplot2$position, y = testplot2$CIL, lty = 'dashed', col = 'blue', lwd = 4)
# lines(x= testplot2$position, y =  testplot2$b_est, type = 'b', col = 'green', pch = 1, lwd = 2)
# lines(x= testplot2$position, y =  testplot2$CIU, lty = 'dashed', col = 'blue', lwd = 4)
# axis(1,                         # Define x-axis manually
#      at = c(0, 50,100, 150, 200, 250),
#      labels = c("1", "50", "100", "150", "200", "averaged 0s'"))
# 
# # Adding a legend
# #legend("bottomright", lty = c(1, 1, 2), legend = c("Truth", "Estimate", "95% CI"),
# #       lwd = 3, col = c("red", "green", "blue"))
# 
# # legend("topright",
# #        # value depending on the windows size
# #        lty = c(1, 1, 2), legend = c("Truth", "Estimate", "95% CI"),
# #        lwd = 3, col = c("red", "green", "blue"),
# #        xpd =FALSE, horiz = TRUE) # You need to specify this graphical parameter to
# # # put the legend outside the plot
# 
# 
# 
# # Back to the default graphical parameters
# #on.exit(par(opar))
# 
