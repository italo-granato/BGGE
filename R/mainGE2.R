#' Genotype x Environment models using linear or gaussian kernel
#'
#' @usage mainGE <- function(Y, X, XF=NULL, W=NULL, method=c("GK", "G-BLUP"), h=NULL, model = c("SM", "MM", "MDs", "MDe", "Cov"), nIter = 150, burnIn = 50, thin = 5, ...)
#'
#' @param Y \code{data.frame} Phenotypic data with three columns. The first column is a \code{factor} for assigned environments,
#' the second column is a \code{factor} for assigned individuals and the third column contains the trait of interest.
#' @param X \code{matrix} Marker matrix with individuals in rows and marker in columns
#' @param XF \code{matrix} Design matrix (\eqn{n \times p}) for fixed effects
#' @param W \code{matrix} Environmental covariance matrix. If \code{W} is provided, it can be used along with \code{MM} and \code{MDs} models. See details
#' @param method Kernel to be used. Methods implemented are the gaussian kernel \code{Gk-EB} and the linear kernel \code{G-BLUP}
#' @param h \code{numeric} Bandwidth parameter to create the gaussian kernel matrix. For \code{method} \code{Gk-EB}, if \code{h} is not provided,
#' then it is computed following a empirical bayesian selection method. See details
#' @param model Specifies the genotype by environment model to be fitted. \code{SM} is the single-environment main genotypic effect model,
#' \code{MM} is the multi-environment main genotypic effect model, \code{MDs} is the multi-environment single variance genotype × environment deviation model,
#' \code{MDe} is the multi-environment, environment-specific variance genotype × environment deviation model
#' @param nIter \code{integer} Number of iterations.
#' @param burnIn \code{integer} Number of iterations to be discarded as burn-in.
#' @param thin \code{integer} Thinin interval.
#' @param \dots additional arguments to be passed.
#' @details
#' The goal is to fit genomic prediction models including GxE interaction. These models can be adjusted through two different kernels.
#' \code{G-BLUP} creates a linear kernel resulted from the cross-product of centered and standarized marker genotypes divide by the number of markers \eqn{p}:
#'     \eqn{G = \frac{XX^T}{p}}
#' If \code{Gk} is choosen, a gaussian kernel is created, resulted from exponential of a genetic distance matrix based on markers scaled by its fifth percentile multiplied by the bandwidth parameter \eqn{h}
#' which has an objective of controlling the rate of decay for correlation between individuals. Thus:
#'  \eqn{ K (x_i, x_{i'}) = exp(-h d_{ii'}^2)}
#' However if the bandwidth parameter is not provided, it need to be estimated from data. The approach currently working is a bayesian method for selecting the bandwidth parameter \eqn{h}
#' through the marginal distribution of \eqn{h}. For more details see Perez-Elizalde et al. (2015).
#' mainGE uses the packages BGLR and MTM to fit the current models:
#' \begin{itemize}
#' \item \code{SM}: is the single-environment main genotypic effect model - The SM model fits the data for each single environment separately.
#' \item \code{MM}: is the multi-environment main genotypic effect model - Multi-environment model considering the main random genetic effects across environments.
#' \item \code{MDs}: is the multi-environment single variance genotype × environment deviation model - This model is an extension of \code{MM} by adding the random interaction effect of the environments
#' with the genetic information of genotypes.
#' \item \code{MDe}: is the multi-environment, environment-specific variance genotype × environment deviation model - This model separates the genetic effects of markers into the main marker effects (across environments) and the specific marker effects (for each environment).
#' \item \code{Cov}: is the multi-environment, variance-covariance environment-specific model -
#'
#' @return
#'
#' @seealso \code{\link[MTM]{MTM}}, \code{\link[BGLR]{BGLR}} and \code{\link{function}}
#'
#'@examples
#'
#'
#'export

mainGE <- function(Y, X, XF=NULL, W=NULL, method=c("GK", "G-BLUP"), h=NULL, model = c("SM", "MM", "MDs", "MDe", "Cov"),
                   nIter = 1000, burnIn = 300, thin = 5, ...) {
  
  names(Y) <- c("environ", "subjects", "value")
  subj <- levels(Y$subjects)
  env <- levels(Y$environ)
  
  #setting kernels
  setK <- getK(Y = Y, X = X, method = method, h = h, model = model)
  K <- setK$K
  y <- as.vector(setK$y)
  
  if (model %in% c("MM", "MDs")) {
    if (hasW) {
      Ze <- model.matrix(~factor(Y[,1])-1)
      Zg <- model.matrix(~factor(Y[,2])-1)
      EC <- Ze %*% (W %*% t(Ze))
      GEC <- ETA$G * ENVc
      K <- c(K, list(EC = list(K = EC), GEC = list(K = GEC)))
    }
  }
  
  if (model == "Cov") {
    tmpY <- matrix(nrow = length(subj), ncol = length(env), data = NA)
    colnames(tmpY) <- environ
    rownames(tmpY) <- subjects
    
    for (i in 1:length(environ)) {
      curEnv <- Y$env == env[i]
      curSub <- match(Y[curEnv, 2], rownames(tmpY))
      tmpY[curSub, i] <- Y[curEnv, 3]
    }
    
    fit <- GEcov(Y = tmpY, K = ETA.tmp$K, nIter = nIter, burnIn = burnIn, thin = thin,...)
  }
  else {
    #' @importFrom BGLR BGLR
    fit <- BGLR(y = Y[,3], ETA = ETA.tmp, nIter = nIter, burnIn = burnIn, thin = thin,...)
  }
  return(fit)
}



getK <- function(Y, X, method = c("GK", "G-BLUP"), h = NULL, model = c("SM", "MM", "MDs", "MDe", "Cov"))
{

  #hasXF <- !is.null(XF)
  subj <- levels(Y[,2])
  env <- levels(Y[,1])
  nEnv <- length(env)
  
  # if(hasXF){
  #   dimXF <- dim(XF)[1]
  #   if ((model == "Cov" &  dimXF != length(subj)) | (model != "Cov" & dimXF != dim(Y)[1]))
  #     stop("Matrix of fixed effects of different dimension")
  # }
  
  if(is.null(rownames(X)))
    stop("Genotype names are missing")
  
  
  if(!all(subj %in% rownames(X)))
    Stop("Not all genotypes presents in phenotypic file are in marker matrix")
  
  X <- X[subj,]
  
  if(model == "SM"){
    if (nEnv > 1)
      stop("Single model choosen, but more than one environment is in the phenotype file")
    
    Zg <- model.matrix( ~ factor(Y[,2L]) - 1)
  }
  else{
    Ze <- model.matrix( ~ factor(Y[,1L]) - 1)
    Zg <- model.matrix( ~ factor(Y[,2L]) - 1)
  }
  
  if (model == "Cov")
    Zg <- model.matrix( ~ factor(subj) - 1)
  
  switch(method,
         'G-BLUP' = {
           # case 'G-BLUP'...
           ker.tmp <- tcrossprod(X) / ncol(X)
           G <- Zg %*% (ker.tmp %*% t(Zg))
         },
         GK = {
           print("Gk")
           # case 'GK'...
           D <- (as.matrix(dist(X))) ^ 2
           naY <- !is.na(Y[, 3])
           Y0 <- Y[naY, 3]
           D0 <- kronecker(matrix(nrow = nEnv, ncol = nEnv, 1), D)
           rownames(D0) <- rep(rownames(D), nEnv)
           D00 <- D0[match(Y[, 2L], rownames(D0)), match(Y[, 2L], rownames(D0))]
           
           if (is.null(h)) {
             sol <- optim(c(1, 1), margh.fun, y = Y0, D = D00, q = quantile(D00, 0.05),
                          method = "L-BFGS-B", lower = c(0.05, 0.05), upper = c(6, 30))
             h <- sol$par[1]
           }
           
           ker.tmp <- exp(-h * D / quantile(D, 0.05))
           G <- Zg %*% (ker.tmp %*% t(Zg))
         },
         {
           stop("Method selected is not available ")
         })
  
  switch(model,
         
         SM = {
           out <- list(K = ker.tmp)
         },
         
         MM = {
           out <- list(G = list(K = G))
         },
         
         MDs = {
           E <- tcrossprod(Ze)
           GE <- G * E
           out <- list(G = list(K = G), GE = list(K = GE))
         },
         
         MDe = {
           ZEE <- matrix(data = 0, nrow = nrow(Ze), ncol = ncol(Ze))
           
           out <- lapply(1:nEnv, function(i){
             ZEE[,i] <- Ze[,i]
             ZEEZ <- ZEE %*% t(Ze)
             K3 <- G * ZEEZ
             return(list(K = K3))
           })
           out <- c(list(list(K = G)), out)
           names(out) <- c("mu", env)
         },
         
         Cov = {
           out <- list(K = ker.tmp)
         }, #DEFAULT CASE
         {
           stop("Model selected is not available ")
         })
  
  
  # ifelse(hasXF, #Test
  #        return(c(list(list(X = XF, model = "FIXED")), out)), #TRUE
         return(list(K = out, y = Y[,3]))#) #FALSE
}


BLMMD <- function(y, K, XF=NULL, me = 1e-10, ite = 10000, burn = 2000, thin = 3, verbose = FALSE) {
  ### PART I  - Conditional distributions functions and eigen descomposition ####
  # Conditional Distribution of tranformed genetic effects b (U'u)
  dcondb <- function(n, media, vari) {
    sd <- sqrt(vari)
    return(rnorm(n, media, sd))
  }
  
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
      Stop("XF is not full-column rank")
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