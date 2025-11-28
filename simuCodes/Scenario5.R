rm(list=ls())
dir.old <-'/home/jiangxl/MultiRFM_scenario5'
setwd(dir.old)

source("helpfunc.R")
source('MSFR_main_R_MSFR_V1.R')

library(dplyr)
library(MultiRFM)
dir.name <- "scenario5"
if(!dir.exists(dir.name)){
  dir.create(dir.name)
}
setwd(dir.name)
# Define functions --------------------------------------------------------
stepwiseSelection <- function(loading, threshold= 1e-3, 
                              method=c("SVR", "var.prop"), greater.than=1, upper.var.prop=0.95){
  
  method <- match.arg(method)
#从给定的选项中选择一个参数值
  qlist <- list()
  
  if(!is.list(loading)){
    d_svdA <- svd(loading)$d
#奇异值分解，打印出奇异值向量
    if(method=='SVR'){
      d_svdA <- d_svdA[d_svdA>threshold]
#保留了所有大于0.001的元素，用于去除噪声或者不重要的特征
      qq <- length(d_svdA)
#计算筛选后的长度
      if(qq>1){
        rsigs <- d_svdA[-qq]/d_svdA[-1]
        qlist$q <-  which.max(rsigs[-c(1:greater.than)]) + greater.than
#移除前greater than个元素，返回最大值索引
      }else{
        qlist$q <- 1
      }
    }else if(method=='var.prop'){
      prop.vec <- cumsum(d_svdA^2)/sum(d_svdA^2)
      qlist$q <- which(prop.vec> upper.var.prop)[1]
    }
    
  }else if(is.list(loading)){
    n_qs <- length(loading)
#算每个特有的因子
    qvec <- rep(NA, n_qs)
    names(qvec) <- paste0("qs", 1:n_qs)
    for(i in 1:n_qs){
      # i <-1 
      d_svdB1 <- svd(loading[[i]])$d
      # 
      if(method=='SVR'){
        d_svdB1 <- d_svdB1[d_svdB1>threshold]
        qq1 <- length(d_svdB1)
        if(qq1>1){
          rsigB1 <- d_svdB1[-qq1]/d_svdB1[-1]
          qvec[i] <-  which.max(rsigB1[-c(1:greater.than)]) + greater.than
        }else{
          qvec[i] <-  1
        }
      }else if(method=='var.prop'){
        prop.vec <- cumsum(d_svdB1^2)/sum(d_svdB1^2)
        qvec[i] <- which(prop.vec> upper.var.prop)[1]
      }
    }
    qlist$qs <- qvec
  }
  
  return(qlist[[1]])
}


getq.from.mat <- function(mat, qmax){
  
  if(all(is.na(mat))) return(c(NA, NA, NA))
  mat[is.na(mat)] <- Inf
  id.min <- which.min(mat)
  hqs <- ceiling(id.min/qmax)
  hq <-  id.min%% qmax
  if(hq==0) hq <- qmax
  return(c(q=hq, qs1=hqs, min.val = mat[id.min]))
}
getq.from.array <- function(arr, qmax){
  
  d3 <- dim(arr)[3]
  mat1 <- matrix(NA, d3, 3)
  for(qs2 in 1:d3){
    # qs2 <- 1
    mat1[qs2, ] <- getq.from.mat(arr[,,qs2], qmax)
  }
  hqs2 <- which.min(mat1[,3])
  
  return(c(q=mat1[hqs2,1], qs1=mat1[hqs2,2], qs2=hqs2))
}


# Select the number of factors --------------------------------------------
library(MultiRFM)
nu.vec <- c(1, 2, 3,  10)
# i <- commandArgs(TRUE) %>% as.integer()
j <- 3
nu <- nu.vec[j]
p <- 100
nvec <- c(150,200);  q <- 3;qs <- c(2,2); S <- length(nvec)
sigma2_eps <- 1
N <- 100
Methods <- c("MultiRFM-SVR", "MSFA-AIC", "MSFA-BIC", "RTS-ER")
nMethods <- length(Methods)
qArray <- array(NA, dim=c(3, nMethods, N))
timeMat <- matrix(NA, N, nMethods)
q_max <- 6; qs_max <- 6

for(i in 1:N){
  # i <- 1
  
  message("i = ", i)
  datList <- gendata_simu_multi(seed=i, nvec=nvec, p=p, q=q, qs=qs, rho=c(6,6), err.type='mvt', sigma2_eps = sigma2_eps, nu=nu)
  # str(datList)
  XList <- datList$Xlist; 
  Lambda  <- matrix(sigma2_eps, p, 2)
  
  #实际的代入计算
  tic <- proc.time()
  res <- MultiRFM(XList, q=q_max, qs= rep(qs_max,2)) 
  hq <- stepwiseSelection(res$A, method='SVR', threshold=1e-10)
  res2 <- MultiRFM(XList, q=hq, qs= rep(qs_max,2)) 
  hq_s <- stepwiseSelection(res2$B, method='SVR')
  toc <- proc.time()
  time.multirfm <- toc[3] - tic[3]
  qArray[,1, i] <-c(hq, hq_s)#
  timeMat[i, 1] <- time.multirfm
  
  
  aicMat1 <- array(NA, dim=c(q_max, qs_max, qs_max))
  bicMat1 <- aicMat1
  maxIter.msfa <- 500
  X_s <- lapply(XList, scale, scale=FALSE)
  hmu <- sapply(XList, colMeans)
  
  tic <- proc.time()
  for(q1 in 1:q_max){
    # q1 <- 1
    message("q1 = ", q1, "/", q_max)
    for(qs1 in 1:qs_max){
      # qs1 <- 1
      for(qs2 in 1:qs_max){
        # qs2 <- 1
        require(MSFA)
        start_value <- start_msfa(X_s =X_s, k = q1, j_s = c(qs1,qs2) )
        mle <-  ecm_msfa(X_s, start_value, trace = TRUE, nIt = maxIter.msfa)
        aicMat1[q1, qs1, qs2] <- mle$AIC
        bicMat1[q1, qs1, qs2] <- mle$BIC
      }
      
    }
  }
  toc <- proc.time()
  time.msfa.log <- toc[3] - tic[3]
  timeMat[i,2] <- time.msfa.log
  timeMat[i,3] <- time.msfa.log+2.5
  qArray[,2,i] <- getq.from.array(aicMat1, qmax = q_max)
  qArray[,3,i] <- getq.from.array(bicMat1, qmax = q_max)
  
  
  # Khat <- SK(Reduce(rbind, X_s))
  tic <- proc.time()
  Khat <- cov(Reduce(rbind, X_s))
  eigvalues  <- eigen(Khat)$values
  nn <- length(eigvalues)
  hq <- which.max(eigvalues[1:(nn-1)] /eigvalues[2:nn] )
  toc <- proc.time()
  qArray[,4, i] <- hq
  timeMat[i,4] <- toc[3] - tic[3]
  
}


print("End!Great Job!")
#saveRDS(qArray, file = "qArrayj=3.rds")
save(qArray, file=paste0("j",j,"N",N,'_qArray.rds'))
