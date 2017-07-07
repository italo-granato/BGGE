#' @title Selection of bandwidth parameter (\eqn{h}) for kernel regression
#' 
#' @description Estimation of bandwidth parameter \eqn{h} of the Gaussian kernel by Bayesian method  
#' 
#'
#' @usage h.fun(Y, D)
#'
#' @param Y \code{data.frame} The first column is a factor denoting environments,
#' the second column is a factor identifying genotypes and the third column contains the trait of interest. Missing values
#' are assigned as NA.
#' @param D \code{matrix}  Genetic distance between individuals.

#' @details
#' The reproducing kernel (RK) function has two components: a genetic distance between individuals based on markers and the bandwidth parameter 
#' \eqn{h} that controls the rate of decay of the covariance between genotypes. Here, \eqn{h} is estimate from data. 
#' The approach is based on the Bayesian method for selecting the bandwidth parameter \eqn{h} through its marginal distribution. For more details see Perez-Elizalde et al. (2015).
#' This function uses all dataset available, however missing data are ignored. If necessary to compute one \eqn{h} for each environment a proper subset of phenotypic data should be used.
#' 
#' @return 
#' it returns the value for \eqn{h}
#' 
#' @references
#' Pérez-Elizalde, S.,Cuevas, J.; Pérez-Rodríguez, P.; Crossa, J. (2015) Selection of The Bandwidth Parameter in a Bayesian Kernel Regression Model for Genomic-Enabled Prediction.
#' J Agr Biol Envir S, 20-4:512-532
#' 
#' 
#' 
#' @examples 
#' # get h for one environment
#' library(BGLR)
#' 
#' data(wheat)
#' X <- scale(wheat.X, scale = TRUE, center = TRUE)
#' D <- (as.matrix(dist(X))) ^ 2
#' rownames(X) <- 1:599
#' pheno_geno <- data.frame(env = gl(n = 1, k = 599), 
#'                GID = gl(n=599, k=1),
#'                value = as.vector(wheat.Y[,1]))
#' h <- h.fun(Y = pheno_geno, D = D)
#' 

#' @export
h.fun <- function(Y, D)
{
  subj <- unique(as.vector(Y[,2]))
  
  if(is.null(rownames(D)))
    stop("Genotype names are missing")
  
  if(!all(subj %in% rownames(D)))
    stop("Not all genotypes presents in phenotypic file are in marker matrix")
  
  nEnv <- length(unique(Y[,1]))
  naY <- !is.na(Y[, 3])
  Y0 <- Y[naY, 3]
  D0 <- kronecker(matrix(nrow = nEnv, ncol = nEnv, 1), D)
  rownames(D0) <- rep(rownames(D), nEnv)
  D00 <- D0[match(Y[, 2L], rownames(D0)), match(Y[, 2L], rownames(D0))]
  D00 <- D00[naY, naY]
  sol <- optim(c(1, 1), margh.fun, y = Y0, D = D00, q = quantile(D00, 0.05),
                 method = "L-BFGS-B", lower = c(0.05, 0.05), upper = c(6, 30))
  h <- sol$par[1]
  
  return(h)
}

margh.fun <- function(theta, y, D, q, nu = 0.0001, Sc = 0.0001, nuh=NULL, Sch=NULL, prior=NULL)
{
  h <- theta[1]
  phi <- theta[2]
  Kh <- exp(-h * D / q)
  eigenKh <- eigen(Kh)
  nr <- length(which(eigenKh$val> 1e-10))
  Uh <- eigenKh$vec[,1:nr]
  Sh <- eigenKh$val[1:nr]
  d <- t(Uh) %*% scale(y, scale = F)
  
  Iden <- -1/2 * sum(log(1 + phi*Sh)) - (nu + nr - 1)/ 2 * log(Sc + sum(d^2 / (1 + phi*Sh)))
  if(!is.null(prior)) lprior <- dgamma(h, nuh, Sch, log = T) else lprior <- 0
  
  Iden <- -(Iden + lprior)
  
  return(Iden)
}