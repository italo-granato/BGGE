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
#' @param verbose Should iteration history be printed on screen? if TRUE or 1 then it is printed,
#' otherwise if another number $n$ is choosen the history is printed every $n$ times. Default is FALSE.  
#' @param me \code{numeric} Trheshold for selecting eigenvalues. Eigenvalues lower than \code{me} are discarded. Default is 1e-10.
#' @details
#' The goal is to fit genomic prediction models including GxE interaction. These models can be adjusted through two different kernels.
#'
#' @examples 
#' # single environment model fit
#' library(BGLR)
#' 
#' data(wheat)
#' X <- scale(wheat.X, scale = TRUE, center = TRUE)
#' K <- list(G = list(K = tcrossprod(X)/ncol(X)))
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
  
  # Conditional fixed effects
  rmvnor <- function(n, media, sigma){
    z <- rnorm(n)
    return(media + crossprod(chol(sigma), z))
  }
  
  
  # Function for eigen descompositions
  eig <- function(K, me) {
    ei <- eigen(K)
    fil <- which(ei$values > me)
    return(list(ei$values[fil], ei$vectors[, fil]))
  }
  
  # verbose part I
  if(as.numeric(verbose) != 0){
    cat("Setting parameters...", "\n", "\n")
  }
  
  ### PART II  Preparation for Gibbs sampler ######
  # Identification of NA's values
  y <- as.numeric(y)
  yNA <- is.na(y)
  whichNa <- which(yNA)
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
    Beta <- solve(crossprod(XF), crossprod(XF, y))
    nBeta <- length(Beta)
    tXX <- solve(crossprod(XF))
  }
  
  # Eigen descomposition for nk Kernels
  nk <- length(K)
  nr <- numeric(nk)
  
  Ei <- list() 
  
  for (i in 1:nk) {
    Ei[[i]] <- list()
    ei <- eig(K[[i]]$K, me)
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
  sigsq.mcmc <- numeric(ite)
  mu.mcmc <- numeric(ite)
  if(!is.null(XF)) betaa.mcmc <- matrix(0, nrow = ite, ncol = nBeta)
  
  sigb.mcmc <- matrix(0, nrow = ite, ncol = nk)
  #u.mcmc <- array(0, dim = c(ite, n, nk))
  u.mcmc <- list()
  u.mcmc[seq(nk)] <- list(matrix(NA_integer_, nrow = ite, ncol = n))
  
  temp <- y - mu
  if(!is.null(XF)){
    temp <- temp - XF %*% beta
  }
  temp <- temp - Reduce('+', u)
  
  i <- 1
  
  ### PART III  Fitted model with training data ####
  
  # Iterations of Gibbs sampler
  while (i <= ite) {
    time.init <- proc.time()[3]
    
    # Conditional of mu
    temp <- temp + mu
    mu <- rnorm(1, mean(temp), sqrt(sigsq/n))
    mu.mcmc[i] <- mu
    temp <- temp - mu
    
    # Conditional of fixed effects
    if(!is.null(XF)){
      temp <- temp + XF %*% Beta
      vari <- sigsq * tXX
      media <- tXX %*% crossprod(XF, temp)
      Beta <- rmvnor(nBeta, media, vari)
      betaa.mcmc[i,] <- Beta
      temp <- temp - XF %*% Beta
    }
    
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
      #u.mcmc[i, , j] <- u[[j]]
      u.mcmc[[j]][i,] <- u[[j]]
    }
    
    # Sampling residual variance 
    res <- temp
    sigsq <- dcondsigsq(res, n, nu, Sce)
    tau <- 1 / sigsq
    sigsq.mcmc[i] <- sigsq
    if(nNa > 0){
      uhat <- Reduce('+', u)
      
      if(!is.null(XF)){
        aux <- XF[yNA,] %*% Beta
      }else{
        aux <- 0
      }
      
      y[whichNa] <- aux + mu + uhat[whichNa] + rnorm(n = nNa, sd = sigsq)
      temp[whichNa] <- y[whichNa] - uhat[whichNa] - aux - mu
    }
    
    if(as.numeric(verbose) != 0 & i %% as.numeric(verbose) == 0){
      time.end <- proc.time()[3]
      cat("Iter: ", i, "time: ", round(time.end - time.init, 3),"\n", "\n")
    }
    i <- i + 1
  }
  
  
  ###### PART IV  Output ######
  #Sampling
  draw <- seq(burn, ite, thin)
  
  # Burning and thinning chains
  sigsq.est <- mean(sigsq.mcmc[draw])
  sigu.est <- apply(sigb.mcmc, MARGIN = 2, FUN = function(x) mean(x[draw]))
  mu.est <- mean(mu.mcmc[draw])

  u.est <- matrix(0, n, nk)
  yHat <- mu.est
  if (!is.null(XF)){
    betaa.est <- colMeans(betaa.mcmc[draw,])
    yHat <- yHat + XF %*% betaa.est
  }
  
  u.est <- sapply(u.mcmc, FUN = function(x) colMeans(x[draw, ]) )
  yHat <- yHat + rowSums(u.est)
  
  result <- list(yHat = yHat, varE = sigsq.est, varU = sigu.est)
  
  return(result)
}
