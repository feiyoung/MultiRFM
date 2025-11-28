# generate man files
# devtools::document()
# R CMD check --as-cran MultiRFM_1.1.0.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("simu_low_dim")
# pkgdown::build_article("ProFASTdlpfc2")


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
#' Fit the high-dimensional multi-study robust factor model
#' @description Fit the high-dimensional multi-study robust factor model which learns latent features and accounts for the heterogeneity among source.
#' @param XList A length-M list, where each component represents a matrix and is the
#' @param q an optional integer, specify the number of study-shared factors; default as 15.
#' @param qs a integer vector with length M, specify the number of study-specifed factors; default as 2.
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization;default as 1.
#' @return return a list including the following components:(1) F, a list composed by the posterior estimation of study-shared factor matrix for each study; (2) H,  a list composed by the posterior estimation of study-specified factor matrix for each study;
#' (3) Sf, a list consisting of the posterior estimation of covariance matrix of study-shared factors for each study; (4) Sh, a list consisting of the posterior estimation of covariance matrix of study-specified factors for each study;
#' (5) A, the loading matrix corresponding to study-shared factors; (6) B, a list composed by the loading matrices corresponding to the study-specified factors;
#' (7) mu,the mean of XList;(8) ELBO: the ELBO value when algorithm stops; (9) ELBO_seq: the sequence of ELBO values.
#' (10) time_use, the elapsed time for model fitting.
#'
#' @details None
#' @export
#' @useDynLib MultiRFM, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#'
#' @examples
#' p <- 100
#' nvec <- c(150,200); qs <- c(2,2)
#' datList <- gendata_simu_multi(seed=1, nvec=nvec, p=p, q=3, qs=qs, rho=c(5,5),
#'         err.type='mvt', sigma2_eps = 1, nu=3)
#' XList <- datList$Xlist;
#' res <- MultiRFM(XList, q=3, qs= qs)
#' str(res)


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


#' Select the number of factors
#' @description Select the number of factors that are shared among studies q and thos that are specific to individual studies(qs).More details are in Section 3.1 of the article.
#' @param XList A length-M list, where each component represents a matrix and is the
#' @param q_max an optional integer, specify the maximum number of study-shared factors; default as 15.
#' @param qs_max an optional integer, specify the maximum number of study-specified factors; default as 4.
#' @param epsELBO  an optional positive value, tolerance of relative variation rate of the evidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization;default as 1.
#' @param threshold the cutoff of the singular values, where the singular values less than this value will be removed.
#' @param method an optional character, contains the methods of "SSVR" and "CUP", where `SSVR` is the sequential singular value ratio method while `CUP` is the criterion based on cumulative proportion of explained variance.
#' @param cup.upper upper limit of the cumulative proportion of explained variance.
#'
#' @return return a list contains the following components:(1) q, the number of shared factors; (2) qs,the number of specified factors.
#'
#' @details None
#' @export
#'
#' @examples
#' p <- 100
#' nvec <- c(150,200); qs <- c(2,2)
#' datList <- gendata_simu_multi(seed=1, nvec=nvec, p=p, q=3, qs=qs, rho=c(5,5),
#'         err.type='mvt', sigma2_eps = 1, nu=3)
#' XList <- datList$Xlist;
#' ## Set maxIter=5 for demonstration while set it to 30 in the formal run.
#' hqlist <- selectFac.MultiRFM(XList, q_max=6, qs_max= rep(4,2), maxIter = 5) #
#' str(hqlist)
#'
#'
#'
selectFac.MultiRFM <- function(XList,  q_max = 15, qs_max = 4, method = c("SSVR", "CUP"),
                               threshold=1e-5, cup.upper=0.95,
                               epsELBO = 1e-05, maxIter = 30, verbose = TRUE, seed = 1){

  S <- length(XList)
  method <- match.arg(method)
  res <- MultiRFM(XList, q=q_max, qs= rep(qs_max,S), epsELBO = epsELBO, maxIter = maxIter, verbose = verbose,
                  seed = seed)
  hq <- stepwiseSelection(res$A, method=method, threshold=threshold, cup.upper=cup.upper)

  res2 <- MultiRFM(XList, q=hq, qs= rep(qs_max,S),  epsELBO = epsELBO, maxIter = maxIter, verbose = verbose,
                   seed = seed)
  hq_s <- stepwiseSelection(res2$B, method=method, threshold=threshold, cup.upper=cup.upper)
  return(list(q=hq, qs = hq_s))
}

stepwiseSelection <- function(loading, threshold= 1e-3,
                              method=c("SSVR", "CUP"), greater.than=1, cup.upper=0.95){

  method <- match.arg(method)
  qlist <- list()

  if(!is.list(loading)){
    d_svdA <- svd(loading)$d
    if(method=='SSVR'){
      d_svdA <- d_svdA[d_svdA>threshold]
      qq <- length(d_svdA)
      if(qq>1){
        rsigs <- d_svdA[-qq]/d_svdA[-1]
        qlist$q <-  which.max(rsigs[-c(1:greater.than)]) + greater.than
      }else{
        qlist$q <- 1
      }
    }else if(method=='CUP'){
      prop.vec <- cumsum(d_svdA^2)/sum(d_svdA^2)
      qlist$q <- which(prop.vec> cup.upper)[1]
    }

  }else if(is.list(loading)){
    n_qs <- length(loading)
    qvec <- rep(NA, n_qs)
    names(qvec) <- paste0("qs", 1:n_qs)
    for(i in 1:n_qs){
      d_svdB1 <- svd(loading[[i]])$d
      if(method=='SSVR'){
        d_svdB1 <- d_svdB1[d_svdB1>threshold]
        qq1 <- length(d_svdB1)
        if(qq1>1){
          rsigB1 <- d_svdB1[-qq1]/d_svdB1[-1]
          qvec[i] <-  which.max(rsigB1[-c(1:greater.than)]) + greater.than
        }else{
          qvec[i] <-  1
        }
      }else if(method=='CUP'){
        prop.vec <- cumsum(d_svdB1^2)/sum(d_svdB1^2)
        qvec[i] <- which(prop.vec> cup.upper)[1]
      }
    }
    qlist$qs <- qvec
  }

  return(qlist[[1]])
}



