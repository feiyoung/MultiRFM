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

approxPCA <- function(X, q){ ## speed the computation for initial values.
  require(irlba)
  n <- nrow(X)
  svdX  <- irlba(A =X, nv = q)
  PCs <- svdX$u * sqrt(n)
  loadings <- svdX$v %*% diag(svdX$d[1:q]) /sqrt(n)
  errMat <- X - PCs %*% t(loadings)
  return(list(PCs = PCs, loadings = loadings, errMat=errMat))
}
normlize <- function(Z){
  nc <- ncol(Z)
  A <- qr(Z)
  A1 <- qr.Q(A)
  A1 <- A1 %*% Diag(sign(A1[1,]))
  return(A1)
}
mat2list <- function(z_int, nvec){

  zList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){

    zList_int[[i]] <- z_int[istart: sum(nvec[1:i]), ]
    istart <- istart + nvec[i]
  }
  return(zList_int)
}
dimrandomfun_norm <- function(q, n){
  matrix(rnorm(prod(q*n)), nrow=n, ncol= q) ## sometimes unform is better, sometimes normal is better!!!
}

MultiRFM <- function (XList,  q = 15, qs = rep(2, length(XList)), epsELBO = 1e-05, maxIter = 30, verbose = TRUE,
                      seed = 1) {

  S <- length(XList); p <- ncol(XList[[1]]); nvec <- sapply(XList, nrow)
  bmu_int <- sapply(XList, colMeans)
  LambdaMat_int<- matrix(1, p, S); nu_int <- 5;
  Xmat <- scale(Reduce(rbind, XList), scale=FALSE);
  fit_pca <- approxPCA(Xmat, q)
  A_int <- fit_pca$loadings
  # SList_f_int <- lapply(nvec, function(x) matrix(1, x, q))
  MuList_f_int <- mat2list(fit_pca$PCs, nvec)

  # SList_h_int <- lapply(1:S, function(s) matrix(1, nvec[s], qs[s]))
  SList_h_int <- SList_f_int <- list()
  for(s in 1:S){
    SList_f_int[[s]] <- array(0, dim=c(q,q, nvec[s]))
    SList_h_int[[s]] <- array(0, dim=c(qs[s],qs[s], nvec[s]))
    for(i in 1:nvec[s]){
      SList_f_int[[s]][,,i] <- diag(rep(1, q))
      SList_h_int[[s]][,,i] <- diag(rep(1, qs[s]))
    }
  }
  set.seed(seed) # random initialization
  MuList_h_int <- lapply(1:S, function(s) matrix(rnorm(nvec[s]* qs[s]), nvec[s], qs[s]))
  BList_int <- lapply(qs, dimrandomfun_norm, n=p)
  BList_int <- lapply(BList_int, normlize)
  tic <- proc.time()
  reslist <- rlfm_cpp(XList, bmu_int, A_int, BList_int, nu_int, LambdaMat_int,
                      SList_f_int, MuList_f_int, SList_h_int, MuList_h_int, epsELBO,
                      maxIter, verbose, loop_ic = FALSE)
  toc <- proc.time()
  time_use <- toc[3] - tic[3]
  reslist$time_use <- time_use
  return(reslist)

}

