
## https://github.com/blhansen/VI-MSFA
rm(list=ls())
dir.old <-'/home/jiangxl/MultiRFM_scenario3'
setwd(dir.old)
source("helpfunc.R")
source('MSFR_main_R_MSFR_V1.R')
library(dplyr)
library(MultiRFM)
dir.name <- "scenario3"
if(!dir.exists(dir.name)){
  dir.create(dir.name)
}
setwd(dir.name)
# Estimation performance of parameters: Fix rhoA=1, change rhoB ------------------------------------
## Fix rhoA=1, change rhoB
nu.vec <- c(2, 3, 20)
# i <- commandArgs(TRUE) %>% as.integer()
j <- 2
nu <- nu.vec[j]
p <- 100
nvec <- c(150,200); d <- 3; q <- 3;qs <- c(2,2); S <- length(nvec)
sigma2_eps <- 1
N <- 100

method_run <- 'All'
methodNames <- c("MultiRFM", "MSFA", "BMSFA", 'MSFR', "MTS", "MSFA-CAVI", "MSFA-SVI")
n_methods <- length(methodNames)
metricList <- list(A_tr=matrix(NA,N, n_methods), B_tr=matrix(NA, N, n_methods),
                   mu_er=matrix(NA, N, n_methods),
                   labmda_er = matrix(NA, N, n_methods),
                   F_tr =matrix(NA,N, n_methods), H_tr =matrix(NA,N, n_methods), 
                   timeMat = matrix(NA, N, n_methods))
for(ii in 1:length(metricList)) colnames(metricList[[ii]]) <- methodNames

for(i in 1:N){
#i <- 1
print("MultiRFM")
message("i = ", i)
datList <- gendata_simu_multi(seed=i, nvec=nvec, p=p, q=q, qs=qs, rho=c(3,3), err.type='mvt', sigma2_eps = sigma2_eps, nu=nu)
# str(datList)
XList <- datList$Xlist; 
Lambda  <- matrix(sigma2_eps, p, 2)

res <- MultiRFM(XList, q=q, qs= qs) 
metricList$timeMat[i,1] <- res$time_use
metricList$A_tr[i,1] <- trace_statistic_fun(res$A, datList$A0)
metricList$B_tr[i,1] <- trace_list_fun(res$B, datList$Blist0)
metricList$F_tr[i,1] <- trace_list_fun(res$F, datList$Flist)
metricList$H_tr[i,1] <- trace_list_fun(res$H, datList$Hlist)
metricList$mu_er[i, 1] <- normvec(datList$mu-res$mu)
metricList$labmda_er[i, 1] <- normvec(Lambda -res$Lambda)



maxIter.msfa <- 10000
## MSFA
print("MSFA")
message("i = ", i)
X_s <- lapply(XList, scale, scale=FALSE)
hmu <- sapply(XList, colMeans)
#remotes::install_github("blhansen/vi-msfa", auth_token = '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')
require(MSFA)
start_value <- start_msfa(X_s =X_s, k = q, j_s = qs)
tic <- proc.time()
mle <-  ecm_msfa(X_s, start_value, trace = TRUE, nIt = maxIter.msfa)
toc <- proc.time()
time.use <- toc[3] - tic[3]
message("Time is :", time.use)
str(mle)
metricList$timeMat[i,2] <- time.use
metricList$A_tr[i,2] <- trace_statistic_fun(mle$Phi, datList$A0)
metricList$B_tr[i,2] <- trace_list_fun(mle$Lambda_s, datList$Blist0)
metricList$mu_er[i, 2] <- normvec(datList$mu-hmu)
metricList$labmda_er[i, 2] <- normvec(Lambda -Reduce(cbind, mle$psi_s))
res_fac_msfa <- estimat.facs(XList, hmu, mle$Phi, mle$Lambda_s)
metricList$F_tr[i,2] <-trace_list_fun(res_fac_msfa$F, datList$Flist)
metricList$H_tr[i,2] <-trace_list_fun(res_fac_msfa$H, datList$Hlist)

# BMSFA
print("BMSFA")
set.seed(1971)
tic <- proc.time()
control.list <- sp_msfa_control(nrun=maxIter.msfa, burn=round(maxIter.msfa/2))
out10_1010 <- sp_msfa(X_s,  k = q,  j_s = qs, trace = TRUE, control = control.list)
toc <- proc.time()
time.use_msfa.bayes <- toc[3] - tic[3]
hA <- apply(out10_1010$Phi, c(1, 2), median)
hBlist <- lapply(out10_1010$Lambda, function(x) apply(x, c(1,2), median))
hLam <- sapply(out10_1010$psi, function(x) apply(x, c(1,2), median))
metricList$timeMat[i,3] <- time.use_msfa.bayes
metricList$A_tr[i,3] <- trace_statistic_fun(hA, datList$A0)
metricList$B_tr[i,3] <- trace_list_fun(hBlist, datList$Blist0)
metricList$mu_er[i, 3] <- normvec(datList$mu-hmu)
metricList$labmda_er[i, 3] <- normvec(Lambda -hLam)
res_fac_msfa.bayes <- estimat.facs(XList, hmu, hA, hBlist)
metricList$F_tr[i,3] <-trace_list_fun(res_fac_msfa.bayes$F, datList$Flist)
metricList$H_tr[i,3] <-trace_list_fun(res_fac_msfa.bayes$H, datList$Hlist)



## RTS
print("RTS")
fit.rts <- RTS.run(X_s, q)
metricList$timeMat[i,5] <- fit.rts$time.use
metricList$A_tr[i,5] <- trace_statistic_fun(fit.rts$Lhat, datList$A0)
metricList$F_tr[i,5] <- trace_list_fun(fit.rts$Flist, datList$Flist)
metricList$mu_er[i, 5] <- normvec(datList$mu-hmu)

library(VIMSFA)
### MSFA-CAVI
print("MSFA-CAVI")
tic <- proc.time()
cavi_est <- cavi_msfa(X_s,  K=q, J_s=qs)
toc <- proc.time()
time_cavi <- toc[3] - tic[3]
hLam <- Reduce(cbind, cavi_est$mean_psi_s)
hF_cavi <- hH_cavi <- list()
for(s in 1:S){
  # s <- 1
  hF_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_f[[s]]))
  hH_cavi[[s]] <- t(Reduce(cbind, cavi_est$mean_l[[s]]))
}

metricList$timeMat[i,6] <- time_cavi
metricList$A_tr[i,6] <- trace_statistic_fun(cavi_est$mean_phi, datList$A0)
metricList$B_tr[i,6] <- trace_list_fun(cavi_est$mean_lambda_s, datList$Blist0)
metricList$mu_er[i,6] <- normvec(datList$mu-hmu)
metricList$labmda_er[i,6] <- normvec(Lambda -hLam)
metricList$F_tr[i,6] <- trace_list_fun(hF_cavi, datList$Flist)
metricList$H_tr[i,6] <- trace_list_fun(hH_cavi, datList$Hlist)

### MSFA-SVI
print("MSFA-SVI")
tic <- proc.time()
svi_est <- svi_msfa(X_s, K=q, J_s=qs)
toc <- proc.time()
time_svi <- toc[3] - tic[3]
hLam <- Reduce(cbind, svi_est$mean_psi_s)
hF_cavi <- hH_cavi <- list()
for(s in 1:S){
  # s <- 1
  hF_cavi[[s]] <- t(Reduce(cbind, svi_est$mean_f[[s]]))
  hH_cavi[[s]] <- t(Reduce(cbind, svi_est$mean_l[[s]]))
}
metricList$timeMat[i,7] <- time_svi
metricList$A_tr[i,7] <- trace_statistic_fun(svi_est$mean_phi, datList$A0)
metricList$B_tr[i,7] <- trace_list_fun(svi_est$mean_lambda_s, datList$Blist0)
metricList$mu_er[i,7] <- normvec(datList$mu-hmu)
metricList$labmda_er[i,7] <- normvec(Lambda -hLam)
metricList$F_tr[i,7] <- trace_list_fun(hF_cavi, datList$Flist)
metricList$H_tr[i,7] <- trace_list_fun(hH_cavi, datList$Hlist)


save(metricList, file=paste0("scenario3",dir.name,"_",method_run,"nsum", sum(nvec),"nu",paste0(nu, collapse = "_"),
                             "rho=c(3,3)","q", q ,"p",p,"N",N, '_metricList.rds'))

}

#metricList$mu_er
#metricList$A_tr
print("rho=c(3,3)")
result_mean=sapply(metricList, colMeans, na.rm=T)
result_sd=sapply(metricList, colSD)
print(result_mean)
print(result_sd)
print("End!Great Job!")
