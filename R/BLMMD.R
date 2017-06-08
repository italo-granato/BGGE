#' Genotype x Environment models using linear or gaussian kernel
#'
#' @usage BLMMD(y, K, XF = NULL, ite = 1000, burn = 200, thin = 3, verbose = FALSE, me = 1e-10)
#'
#' @param y vector of phenotypic data 
#' @param K Kernels with effects to be fitted
#' @param XF \code{matrix} Design matrix (\eqn{n \times p}) for fixed effects
#' @param ite \code{integer} Number of iterations.
#' @param burn \code{integer} Number of iterations to be discarded as burn-in.
#' @param thin \code{integer} Thinin interval.
#' @param verbose \code{logical} Should report be printed on screen?
#' @param me \code{numeric} tolerance for zero. Default is 1e-10
#' @details
#' The goal is to fit genomic prediction models including GxE interaction. These models can be adjusted through two different kernels.
#'
#' @examples 
#' # single environment model fit
#' library(BGLR)
#' 
#' data(wheat)
#' X <- scale(wheat.X, scale = TRUE, center = TRUE)
#' K <- list(tcrossprod(X)/ncol(X))
#' y <- as.vector(wheat.Y[,1])
#' 
#' fit <- BLMMD(y=y, K=K)
#'
#' @seealso 
#' \code{\link[BGLR]{BGLR}}
#' 
#' 
#' @export
BLMMD <- function(y, K, XF = NULL, ite = 1000, burn = 200, thin = 3, verbose = FALSE, me = 1e-10) {
  ### PART I  - Conditional distributions functions and eigen descomposition ####
  # Conditional Distribution of tranformed genetic effects b (U'u)
  #' @import stats
  dcondb <- function(n, media, vari) {
    sd <- sqrt(vari)
    return(rnorm(n, media, sd))
  }
  
  #' @importFrom stats rgamma
  # Conditional distribution of compound variance of genetic effects (sigb)
  dcondsigb <- function(b, deltav, n, nu, Sc) {
    z <- sum(b ^ 2 * deltav)
    return(1 / rgamma(1, (n + nu) / 2, (z + nu * Sc) / 2))
  } 
  
  # Conditional distribution of scale hyperparameter chi-square (Sc) of 
  # compound variance of genetic effects. 
  dcondSc <- function(nu, sigb) {
    return(rgamma(1, (nu + 1) / 2, (nu / (2 * sigb))))
  }
  
  # Conditional distribution of residual compound variance 
  dcondsigsq <- function(Aux, n, nu, Sce) {
    return(1 / rgamma(1, (n + nu) / 2, crossprod(Aux) / 2 + Sce / 2))
  }
  
  
  # Function for eigen descompositions
  eig <- function(K, me) {
    ei <- eigen(K)
    fil <- which(ei$values > me)
    return(list(ei$values[fil], ei$vectors[, fil]))
  }
  ### PART II  Preparation for Gibbs sampler ######
  
  # Identification of NA's values
  whichNa <- which(is.na(y))
  nNa <- length(whichNa)
  mu <- mean(y, na.rm = TRUE)
  n <- length(y)
  
  if(nNa > 0){ 
    y[whichNa] <- mu
  }
  
  # initial values of fixed effects 
  if(!is.null(XF)) {
    rankXF <- qr(XF)$rank
    if (rankXF < ncol(XF) ){
      stop("XF is not full-column rank")
    }
    beta <- solve(crossprod(XF), crossprod(XF, y))
    if(nNa > 0){
      XFt <- XF[whichNa,]
    }
  }
  
  # Eigen descomposition for nk Kernels
  nk <- length(K)
  nr <- numeric(nk)
  
  Ei<-list() 
  
  for (i in 1:nk) {
    Ei[[i]] <- list()
    ei <- eig(K[[i]], me)
    Ei[[i]]$s <- ei[[1]]
    Ei[[i]]$U <- ei[[2]]
    Ei[[i]]$tU <- t(ei[[2]])
    nr[i] <- length(Ei[[i]]$s)
  }
  
  
  # Initial values for Monte Carlo Markov Chains (MCMC)
  nu <- 3
  Sce <- (nu + 2) * var(y) / 2
  tau <- 0.01
  u <- list()
  for (j in 1:nk) u[[j]] <- rnorm(n, 0, 1 / (2 * n))
  sigsq <- var(y)
  Sc <- rep(0.01, nk)
  sigb <- rep(0.2, nk)
  sigsq.mcmc <- rep(0, ite)
  sigb.mcmc <- matrix(0, nrow = ite, ncol = nk)
  u.mcmc <- array(0, dim = c(ite, n, nk))
  temp <- y
  if(!is.null(XF)){
    temp <- temp - XF %*% beta
  }
  for (k in 1:nk)  temp <- temp - u[[k]]
  
  i <- 1
  
  ### PART III  Fitted model with training data ####
  
  # Iterations of Gibbs sampler
  while (i <= ite) {
    time.init <- proc.time()[3]
    
    # Conditionals x Kernel
    for (j in 1:nk) {
      # Sampling genetic effects
      temp <- temp + u[[j]]
      d <- crossprod(Ei[[j]]$U, temp)
      s <- Ei[[j]]$s
      deltav <- 1 / s
      lambda <- sigb[j]
      vari <- s * lambda / (1 + s * lambda * tau)
      media <- tau * vari * d
      b <- (dcondb(nr[j], media, vari))
      u[[j]] <-crossprod(Ei[[j]]$tU ,b)
      temp <- temp - u[[j]]
      
      # Sampling scale hyperparameters and variance of genetic effects
      Sc[j] <- dcondSc(nu, sigb[[j]])
      sigb[j] <- dcondsigb(b, deltav, nr[j], nu, Sc[j])
      sigb.mcmc[i, j] <- sigb[j]
      u.mcmc[i, , j] <- u[[j]] 
    }
    
    # Sampling residual variance 
    res <- temp
    sigsq <- dcondsigsq(res, n, nu, Sce)
    tau <- 1 / sigsq
    sigsq.mcmc[i] <- sigsq
    if(nNa > 0){
      if(!is.null(XF)){
        XFna <- XFt %*% beta
      }else{
        XFna <- 0
      }
      uhat <- Reduce('+', u)
      y[whichNa] <- XFna + uhat[whichNa] + rnorm(n = nNa, sd = sigsq)
      temp[whichNa] <- y[whichNa] - uhat[whichNa] - XFna
    }
    
    if(verbose){
      time.end <- proc.time()[3]
      cat("iter: ", i, "time: ", round(time.end-time.init, 3),"\n", "\n")
    }
    i <- i + 1
  }
  
  
  ###### PART IV  Output ######
  #Sampling
  draw <- seq(burn, ite, thin)
  
  # Burning and thinning chains
  sigsq.est <- mean(sigsq.mcmc[draw])
  if(class(sigb.mcmc[draw,]) != "matrix"){
    sigu.est <- mean(sigb.mcmc[draw, ])
  }else{
    sigu.est <- apply(sigb.mcmc[draw, ], 2, mean)
  }
  
  u.est <- matrix(0, n, nk)
  yHat <- mu
  if (!is.null(XF)){
    yHat <- XF %*% beta
  }
  
  for (i in 1:nk) {
    u.est[, i] <- apply(u.mcmc[draw, , i], 2, mean)
    yHat <- yHat + u.est[, i]
  }
  
  result <- list(yHat = yHat, varE = sigsq.est, varU = sigu.est)
  
  return(result)
}
