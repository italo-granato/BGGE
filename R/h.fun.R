#' Genotype x Environment models using linear or gaussian kernel
#'
#' @usage h.fun(Y, X)
#'
#' @param Y \code{data.frame} Phenotypic data with three columns. The first column is a \code{factor} for assigned environments,
#' the second column is a \code{factor} for assigned individuals and the third column contains the trait of interest.
#' @param X \code{matrix} Marker matrix with individuals in rows and marker in columns

#' @details
#' The goal is to estimate the bandwith parameter from data. The approach used is a bayesian method for selecting the bandwidth parameter \eqn{h}
#' through the marginal distribution of \eqn{h}. For more details see Perez-Elizalde et al. (2015).

#' export
h.fun <- function(Y, X)
{
  nEnv <- length(unique(Y[,1]))
  D <- (as.matrix(dist(X))) ^ 2
  naY <- !is.na(Y[, 3])
  Y0 <- Y[naY, 3]
  D0 <- kronecker(matrix(nrow = nEnv, ncol = nEnv, 1), D)
  rownames(D0) <- rep(rownames(D), nEnv)
  D00 <- D0[match(Y[, 2L], rownames(D0)), match(Y[, 2L], rownames(D0))]
  sol <- optim(c(1, 1), margh.fun, y = Y0, D = D00, q = quantile(D00, 0.05),
                 method = "L-BFGS-B", lower = c(0.05, 0.05), upper = c(6, 30))
  h <- sol$par[1]
  
  return(h)
}

margh.fun <- function(theta, y, D, q, nu=0.0001, Sc=0.0001, nuh=NULL, Sch=NULL, prior=NULL)
{
  h <- theta[1]
  phi <- theta[2]
  Kh <- exp(-h*D/q)
  eigenKh <- eigen(Kh)
  nr <- length(which(eigenKh$val> 1e-10))
  Uh <- eigenKh$vec[,1:nr]
  Sh <- eigenKh$val[1:nr]
  d <- t(Uh) %*% scale(y, scale=F)
  
  Iden <- -1/2*sum(log(1+phi*Sh)) - (nu+nr-1)/2*log(Sc+sum(d^2/(1+phi*Sh)))
  if(!is.null(prior)) lprior <- dgamma(h,nuh,Sch,log=T) else lprior <- 0
  
  Iden <- -(Iden+lprior)
  
  return(Iden)
}