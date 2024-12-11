


# Our method --------------------------------------------------------------
stepwiseSelection <- function(loading, threshold= 1e-3, 
                              method=c("SVR", "var.prop"), greater.than=1,
                              upper.var.prop=0.95, return.svd=FALSE){
  
  method <- match.arg(method)
  qlist <- list()
  
  if(!is.list(loading)){
    d_svdA <- svd(loading)$d
    if(method=='SVR'){
      d_svdA <- d_svdA[d_svdA>threshold]
      qq <- length(d_svdA)
      if(qq>1){
        rsigs <- d_svdA[-qq]/d_svdA[-1]
        qlist$q <-  which.max(rsigs[-c(1:greater.than)]) + greater.than
        if(return.svd){
          qlist$d_svd <- d_svdA
        }
      }else{
        qlist$q <- 1
      }
    }else if(method=='var.prop'){
      prop.vec <- cumsum(d_svdA^2)/sum(d_svdA^2)
      qlist$q <- which(prop.vec> upper.var.prop)[1]
    }
    
  }else if(is.list(loading)){
    n_qs <- length(loading)
    qvec <- rep(NA, n_qs)
    names(qvec) <- paste0("qs", 1:n_qs)
    d_svdlist <- list()
    for(i in 1:n_qs){
      # i <-1 
      d_svdB1 <- svd(loading[[i]])$d
      d_svdlist[[i]] <- d_svdB1
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
    if(return.svd){
      qlist$d_svd <- d_svdlist
    }
  }
  
  return(qlist)
}

# Generat data ------------------------------------------------------------

gendata_simu_multi <-function (seed = 1, nvec = c(100,300), p = 50, q = 3,
                               qs= rep(2, length(nvec)), err.type=c("gaussian", 'mvt', "exp", "t", "mixnorm", "pareto"),
                               rho = c(1,1), sigma2_eps=0.1, nu=1){
  # seed = 1; nvec = c(100,300); p = 50; q = 3
  # qs= rep(2, length(nvec))
  # rho = 1; sigma2_eps=0.1;
  
  
  if(length(nvec)<2) stop("nvec must have at least two elements!")
  err.type <- match.arg(err.type)
  S <- length(nvec)
  require(MASS)
  factor_term_A <- rho[1]
  factor_term_B <- rho[2]
  set.seed(1) 
  Blist <- list()
  set.seed(1)
  bmu0 <- matrix(rnorm(p*S), p, S) 
  Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
  A <- qr(Ztmp)
  A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
  A1 <- A1 %*% Diag(sign(A1[1, ]))*factor_term_A ## Fixed B0 and mu0 for each repeat.
  A0 <- A1[,1:q]; B1 <- A1[,(q+1):(q+qs[1])]
  t(A0) %*% A0; t(B1) %*% B1
  Blist[[1]] <- B1
  for(r in 2:S){
    set.seed(r)
    Ztmp <- matrix(rnorm(p * (qs[2])), p, (qs[2]))
    A <- qr(Ztmp)
    A1 <- qr.Q(A) %*% Diag(seq(qs[2], 1, length=qs[1]))  * factor_term_B
    t(A1) %*% A1
    Blist[[r]]   <- A1 %*% Diag(sign(A1[1, ])) 
  }
  
  
  
  set.seed(seed)
  Xlist <- list()
  Flist <- list()
  Hlist <- list()
  for(s in 1:S){
    
    n <- nvec[s]
    if(err.type == 'gaussian'){
      epsi <- MASS::mvrnorm(n, mu=rep(0, p), Sigma = diag(rep(p,1))* sigma2_eps)
    }else if(err.type == 'mvt'){
      epsi <- mvtnorm::rmvt(n, sigma =sigma2_eps* diag(p), df=nu)
      #matrix(rt(n*p, df=nu), n, p) * 
    }else if(err.type == 'exp'){
      epsi <- (matrix(rexp(n*p), n, p)-1) * sqrt(sigma2_eps)
      #matrix(rt(n*p, df=nu), n, p) * 
    }else if(err.type == 't'){
      epsi <- matrix(rt(n*p, df=nu), n, p)  * sqrt(sigma2_eps)
    }else if(err.type == 'mixnorm'){
      library(mixtools)
      
      lambda <- rep(1, 2)/2
      mu <- c(-5, 5)
      sigma <- rep(1, 2)
      epsi <- matrix(rnormmix(n*p, lambda, mu, sigma), n, p)  * sqrt(sigma2_eps)
    }else if(err.type== 'pareto'){
      library(LaplacesDemon)
      alpha=2;n <- 1000
      x <- rpareto(n*p, alpha) - alpha/(alpha-1)
      epsi <- matrix(x, n, p)  * sqrt(sigma2_eps)
    }
    
    FH <- mvrnorm(n, mu = rep(0, q+qs[s]), Diag(rep(1,q+qs[s])))
    Flist[[s]] <- FH[,1:q]
    Hlist[[s]] <- FH[,(q+1):(q+qs[s])]
    AB1 <- cbind(A0, Blist[[s]])
    Xlist[[s]]  <- matrix(bmu0[,s], nrow=n, ncol=p, byrow = TRUE) + FH %*% t(AB1) + epsi
  }
  
  return(list(Xlist = Xlist, mu0=bmu0, A0 = A0, Blist0 = Blist,
              Flist = Flist, Hlist = Hlist,  q=q, qs=qs))
}



Diag<-function (vec){
  q <- length(vec)
  if (q > 1) {
    y <- diag(vec)
  }
  else {
    y <- matrix(vec, 1, 1)
  }
  return(y)
}

# approxPCA <- function(X, q){ ## speed the computation for initial values.
#   require(irlba) 
#   n <- nrow(X)
#   svdX  <- irlba(A =X, nv = q)
#   PCs <- svdX$u * sqrt(n)
#   loadings <- svdX$v %*% diag(svdX$d[1:q]) /sqrt(n)
#   errMat <- X - PCs %*% t(loadings)
#   return(list(PCs = PCs, loadings = loadings, errMat=errMat))
# }
# normlize <- function(Z){
#   nc <- ncol(Z)
#   A <- qr(Z)
#   A1 <- qr.Q(A)
#   A1 <- A1 %*% Diag(sign(A1[1,])) 
#   return(A1)
# }
# mat2list <- function(z_int, nvec){
#   
#   zList_int <- list()
#   istart <- 1
#   for(i in 1:length(nvec)){
#     
#     zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
#     istart <- istart + nvec[i]
#   }
#   return(zList_int)
# }
# dimrandomfun_norm <- function(q, n){
#   matrix(rnorm(prod(q*n)), nrow=n, ncol= q) ## sometimes unform is better, sometimes normal is better!!!
# }

# MultiRFM <- function (XList,  q = 15, qs = rep(2, length(XList)), epsELBO = 1e-05, maxIter = 30, verbose = TRUE, 
#                   seed = 1) {
#   
#   S <- length(XList); p <- ncol(XList[[1]]); nvec <- sapply(XList, nrow)
#   bmu_int <- sapply(XList, colMeans)
#   LambdaMat_int<- matrix(1, p, S); nu_int <- 5;
#   Xmat <- scale(Reduce(rbind, XList), scale=FALSE);
#   fit_pca <- approxPCA(Xmat, q)
#   A_int <- fit_pca$loadings
#   # SList_f_int <- lapply(nvec, function(x) matrix(1, x, q))
#   MuList_f_int <- mat2list(fit_pca$PCs, nvec)
#   MuList_h_int <- lapply(1:S, function(s) matrix(0, nvec[s], qs[s]))
#   # SList_h_int <- lapply(1:S, function(s) matrix(1, nvec[s], qs[s]))
#   SList_h_int <- SList_f_int <- list()
#   for(s in 1:S){
#     SList_f_int[[s]] <- array(0, dim=c(q,q, nvec[s]))
#     SList_h_int[[s]] <- array(0, dim=c(qs[s],qs[s], nvec[s]))
#     for(i in 1:nvec[s]){
#       SList_f_int[[s]][,,i] <- diag(rep(1, q))
#       SList_h_int[[s]][,,i] <- diag(rep(1, qs[s]))
#     }
#   }
#   set.seed(seed) # random initialization
#   BList_int <- lapply(qs, dimrandomfun_norm, n=p)
#   BList_int <- lapply(BList_int, normlize)
#   tic <- proc.time()
#   reslist <- rlfm_cpp(XList, bmu_int, A_int, BList_int, nu_int, LambdaMat_int, 
#                       SList_f_int, MuList_f_int, SList_h_int, MuList_h_int, epsELBO, 
#                       maxIter, verbose, loop_ic = FALSE) 
#   toc <- proc.time()
#   time_use <- toc[3] - tic[3]
#   reslist$time_use <- time_use
#   return(reslist)
#   
# }

# Compared methods --------------------------------------------------------
mat2list <-function(z_int, nvec){
  
  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}

estimat.facs <- function(XList, hmu, hA, hB){
  FList <- list()
  HList <- list()
  p <- ncol(XList[[1]])
  S <- length(XList)
  for(s in 1:S){
    ns <- nrow(XList[[s]])
    Xms <- (XList[[s]]-matrix(hmu[,1], nrow=ns, ncol=p)) 
    FList[[s]] <- Xms %*% hA %*% qr.solve(t(hA)%*% hA)
    HList[[s]] <- (Xms - FList[[s]]%*% t(hA)) %*% hB[[s]] %*% qr.solve(t(hB[[s]])%*% hB[[s]])
  }
  return(list(F=FList, H=HList))
}

MSFR.run <- function(XList, ZList, q, qs, maxIter=1e4, load.source=FALSE, dir.source=NULL){
  
  # require(MSFA)
  require(psych)
  if(!load.source){
    source(paste0(dir.source, "MSFR_main_R_MSFR_V1.R"))
  }
  #fa <- psych::fa
  B_s <- ZList
 
  X_s <- XList #
  t1<- proc.time()
  test <- start_msfa1(X_s, B_s, 5, k=q, j_s=qs, constraint = "block_lower2", method = "adhoc")
  # EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE, nIt=maxIter, constraint = "block_lower1")
  EM_beta <- ecm_msfa1(X_s, B_s, start=test, trace = FALSE, nIt=maxIter)
  t2<- proc.time()
  EM_beta$Flist <- lapply(EM_beta$E_f, t)
  EM_beta$Hlist <- lapply(EM_beta$E_l, t)
  EM_beta$time.use <- t2[3] - t1[3]
  return(EM_beta)
}


#----------------------------------------------------------------------------------------
#  Spatial Multivariate Kendall' tau Matrix
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#  Output:
#          TK ------ Sample Spatial Multivariate Kendall' tau Matrix
#----------------------------------------------------------------------------------------

SK<-function(X){
  p <- ncol(X)
  n <- nrow(X)
  TK <- matrix(0,p,p)
  for(i in 1:(n-2)){
    TT <- matrix(rep(X[i,],n-i),n-i,p,byrow = TRUE)-X[(i+1):n,]
    TT <- t(diag(1/diag(TT%*%t(TT)))%*%TT)%*%TT
    TK <- TK+TT
  }
  TT <- X[n-1,]-X[n,]
  TK <- TK+TT%*%t(TT)/sum(TT^2)
  TK <- 2/(n*(n-1))*TK
  return(TK)
}

#----------------------------------------------------------------------------------------
#  RTS Method
#  Input:
#           x ------ n x p  data matrix (n: time dimension, p: cross-section)
#           r ------ number of factors
#  Output:
#        Fhat ------ n x r  the estimated factor scores 
#        Lhat ------ p x r  the estimated factor loadings 
#----------------------------------------------------------------------------------------

RTS <- function(X,r){
  p <- ncol(X)
  n <- nrow(X)
  Khat <- SK(X)
  Lhat <- sqrt(p)*as.matrix(eigen(Khat)$vectors[,1:r]) 
  Fhat <- matrix(0,n,r)
  for (i in 1:n){
    Fhat[i,]=lm(X[i,]~Lhat-1)$coefficients
  }
  return(list(Fhat=Fhat,Lhat=Lhat))
}

### He Yong, 2022, JBES, Robust latent factor model without moment constraint.
RTS.run <- function(XList, q){
  nvec <- sapply(XList, nrow)
  X <- Reduce(rbind, XList)
  tic <- proc.time()
  fit <- RTS(X, r=q)
  toc <- proc.time()
  Flist <- mat2list(fit$Fhat, nvec)
  fit$Flist <- Flist
  fit$time.use <- toc[3] - tic[3]
  return(fit)
}

# Metrics -----------------------------------------------------------------
# mean.Fnorm <- function(x) sum(x^2)/ length(x)
normvec <- function(x) sqrt(sum(x^2)/ length(x))
colSD <- function(X) apply(X, 2, sd, na.rm=TRUE)
trace_statistic_fun <- function(H, H0){
  
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% qr.solve(t(H) %*% H) %*% t(H) %*% H0
  
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
  
}
trace_list_fun <- function(Hlist, H0list){
  trvec <- rep(NA, length(Hlist))
  for(i in seq_along(trvec)){
    trvec[i] <- trace_statistic_fun(Hlist[[i]], H0list[[i]])
  }
  return(mean(trvec))
}

