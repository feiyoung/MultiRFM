
# Generate simulated data -------------------------------------------------
#' Generate Simulated Multi-Study Factor Analysis Data
#'
#' Generate simulated data for multi-study factor analysis under
#' different error distributions. The data follows a factor model with common factors (shared across studies)
#' and study-specific factors (unique to each study), plus noise.
#'
#' @param seed Integer, default = 1. Random seed for reproducibility of simulated data.
#' @param nvec Numeric vector (length >= 2). Sample sizes of each study (e.g., `c(150, 200)` for 2 studies with 150 and 200 samples).
#' @param p Integer, default = 50. Number of variables (features) in the data.
#' @param q Integer, default = 3. Number of common factors (shared across all studies).
#' @param qs Numeric vector with length equal to `length(nvec)`, default = `rep(2, length(nvec))`.
#'   Number of study-specific factors for each study (e.g., `c(2,2)` for 2 studies each with 2 specific factors).
#' @param err.type Character, default = "gaussian". Error distribution type, one of:
#'   - "gaussian": Gaussian (normal) distribution;
#'
#'   - "mvt": Multivariate t-distribution;
#'
#'   - "exp": Exponential distribution (centered to mean 0);
#'
#'   - "t": Univariate t-distribution (independent across variables);
#'
#'   - "mixnorm": Mixture of two normal distributions;
#'
#'   - "pareto": Pareto distribution (centered to mean 0).
#' @param rho Numeric vector of length 2, default = `c(1,1)`. Scaling factors for:
#'   - `rho1`: Common factor loadings (matrix `A0`);
#'   - `rho2`: Study-specific factor loadings (matrix list `Blist0`).
#' @param sigma2_eps Numeric, default = 0.1. Variance of the error term (controls noise level).
#' @param nu Integer, default = 1. Degrees of freedom for t-distribution ("mvt" or "t" `err.type`).
#'   Ignored for other error distributions.
#'
#' @return A list containing the simulated data and true parameter values (for model evaluation):
#' \itemize{
#'   \item{\code{Xlist}: List of matrices. Each element is a data matrix (ns × p) for study s,
#'         where ns = `nvec[s]` (sample size of study s), p = number of variables.}
#'   \item{\code{mu0}: Matrix (p × S). True mean vector for each variable (row) in each study (column),
#'         where S = `length(nvec)` (number of studies).}
#'   \item{\code{A0}: Matrix (p × q). True common factor loadings (shared across all studies) —
#'         constructed as the first q columns of an orthogonal matrix (`A1`) generated internally.
#'         This is the "ground truth" that modeling functions (e.g., MultiRFM) aim to estimate.}
#'   \item{\code{Blist0}: List of matrices. Each element is a true study-specific factor loadings matrix (p × qs[s])
#'         for study s. Constructed from orthogonal matrices (similar to `A0`) and scaled by `rho[2]`.
#'         Another "ground truth" for model evaluation.}
#'   \item{\code{Flist}: List of matrices. Each element is a true common factor score matrix (ns × q) for study s,
#'         generated from a standard normal distribution. These are the latent common factor values used to generate `Xlist`.}
#'   \item{\code{Hlist}: List of matrices. Each element is a true study-specific factor score matrix (ns × qs[s])
#'         for study s, generated from a standard normal distribution. Latent specific factor values used to generate `Xlist`.}
#'   \item{\code{q}: Integer. Number of common factors used for data generation (same as input `q`, for reference).}
#'   \item{\code{qs}: Numeric vector. Number of study-specific factors used for data generation (same as input `qs`, for reference).}
#' }
#'
#' @details
#' The simulated data follows the multi-study factor model:
#'
#' Xs = mu0s + Fs x A0 + Hs x B0s + epsilons
#'
#' True parameters (`A0`, `Blist0`, `mu0`) are generated with orthogonal constraints to ensure identifiability.
#'
#' @examples
#' # Example 1: Gaussian error (2 studies, 100/200 samples, 50 variables)
#' set.seed(123)
#' sim_data <- gendata_simu_multi(
#'   seed = 123,
#'   nvec = c(100, 200),
#'   p = 50,
#'   q = 3,          # 3 common factors
#'   qs = c(2, 2),   # 2 specific factors per study
#'   err.type = "gaussian",
#'   rho = c(1, 1),
#'   sigma2_eps = 0.1
#' )
#' str(sim_data)  # Check structure of simulated data
#'

#' # Extract true parameters for model evaluation
#' true_A <- sim_data$A0        # True common loadings
#' true_B1 <- sim_data$Blist0[[1]]  # True specific loadings (study 1)
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @importFrom mixtools rnormmix
#' @importFrom LaplacesDemon rpareto
#' @importFrom stats rnorm rt rexp
#'
#' @author Wei Liu
#' @export


gendata_simu_multi <-function (seed = 1, nvec = c(100,300), p = 50, q = 3,
                               qs= rep(2, length(nvec)), err.type=c("gaussian", 'mvt', "exp", "t", "mixnorm", "pareto"),
                               rho = c(1,1), sigma2_eps=0.1, nu=1){
  # seed = 1; nvec = c(100,300); p = 50; q = 3
  # qs= rep(2, length(nvec))
  # rho = 1; sigma2_eps=0.1;


  if(length(nvec)<2) stop("nvec must have at least two elements!")
  err.type <- match.arg(err.type)
  S <- length(nvec)
  # require(MASS)
  factor_term_A <- rho[1]
  factor_term_B <- rho[2]
  Blist <- list()
  bmu0 <- matrix(rnorm(p*S), p, S)
  Ztmp <- matrix(rnorm(p * (q+qs[1])), p, (q+qs[1]))
  A <- qr(Ztmp)
  A1 <- qr.Q(A) %*% Diag(seq(q+qs[1], 1, length=q+qs[1]))
  A1 <- A1 %*% Diag(sign(A1[1, ]))*factor_term_A ## Fixed B0 and mu0 for each repeat.
  A0 <- A1[,1:q]; B1 <- A1[,(q+1):(q+qs[1])]
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


      lambda <- rep(1, 2)/2
      mu <- c(-5, 5)
      sigma <- rep(1, 2)
      epsi <- matrix(rnormmix(n*p, lambda, mu, sigma), n, p)  * sqrt(sigma2_eps)
    }else if(err.type== 'pareto'){
      #library(LaplacesDemon)
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

