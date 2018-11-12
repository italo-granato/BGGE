#' Genotype x Environment models using regression kernel
#' 
#' BGGE function fits Bayesian regression for continuous observations through regression kernels
#'
#' @usage BGGE(y, K, XF = NULL, ne, ite = 1000, burn = 200, thin = 3, verbose = TRUE, 
#'             tol = 1e-10, R2 = 0.5)
#'
#' @param y Vector of data. Should be numeric and NAs are allowed.
#' @param K list A two-level list Specify the regression kernels (co-variance matrix). The former is the \code{Kernel},
#' where is included the regression kernels. The later is the \code{Type}, specifying if the matrix is either \code{D} Dense or
#' \code{BD} Block Diagonal. A number of regression kernels or random effects to be fitted are specified in this list.  
#' @param XF matrix Design matrix (\eqn{n \times p}) for fixed effects
#' @param ne vector Number of genotypes by environment. 
#' @param ite numeric Number of iterations.
#' @param burn numeric Number of iterations to be discarded as burn-in.
#' @param thin numeric Thinin interval.
#' @param verbose Should iteration history be printed on console? If TRUE or 1 then it is printed,
#' otherwise, if another number $n$ is choosen the history is printed every $n$ times. The default is \code{FALSE}  
#' @param tol a numeric tolerance level. Eigenvalues lower than \code{tol} are discarded. Default is \code{1e-10}.
#' @param R2 the proportion of variance expected to be explained by the regression.
#' 
#' @details
#' The goal is to fit genomic prediction models for continuous outcomes through Gibbs sampler. BGGE uses a proposal for dimension reduction
#' through an orthogonal transformation of observed data (y) as well as differential shrinkage because of the prior variance assigned 
#' to regression parameters. Further details on this approach can be found in Cuevas et al. (2014).
#' The primaty genetic model is
#' \deqn{y = g + e}
#' where \eqn{y} is the response, \eqn{g} is the unknown random effect and \eqn{e} is the residual effect.
#' You can specify a number of random effects \eqn{g}, as many as desired, through a list of regression kernels related to each random effect in the
#' argument \code{K}.
#' The structure of \code{K} is a two level list, where the first element on the second level is the Kernel and the second element is a definition of
#' type of matrix. There are two definitions, either matrix is \code{D} (dense) or \code{BD} (Block Diagonal). As we make the spectral decomposition 
#' on the kernels, for block diagonal matrices, we take advantage of its structure and make decomposition on the submatrices instead of one big matrix.
#' For example, the regression kernels should be an structure like K = list(list(Kernel = G, Type = "D"), list(Kernel = G, Type = "BD")). 
#' The definition of one matrix as a block diagonal must be followed by the number of subjects in each submatrix in the block diagonal,
#' present in the \code{ne}, which allows sub matrices to be drawn. Some genotype by environment models has the block diagonal matrix type or similar.
#' The genotype x environment deviation matrix in MDs model (Sousa et al., 2017) has the structure of block diagonal. 
#' Also, the matrices for environment-specific variance in MDe models (Sousa et al., 2017) if summed, can form a structure of block diagonal, 
#' where is possible to extract sub matrices for each environment. In the case of all kernel be of the dense type, \code{ne} is ignored. 
#' 
#' @return
#'  A list with estimated posterior means of residual and genetic variance component for each term in the linear model and the genetic value predicted. Also the 
#'  values along with the chains are released. 
#' 
#' 
#' @examples 
#' # multi-environment main genotypic model
#' library(BGLR)
#' data(wheat)
#' X<-wheat.X[1:200,1:600]  # Subset of 200 subjects and 600 markers
#' rownames(X) <- 1:200
#' Y<-wheat.Y[1:200,]
#' A<-wheat.A[1:200,1:200] # Pedigree
#' 
#' GB<-tcrossprod(X)/ncol(X)
#' K<-list(G = list(Kernel = GB, Type = "D"))
#' y<-Y[,1]
#' fit<-BGGE(y = y,K = K, ne = length(y), ite = 300, burn = 100, thin = 2)
#' 
#' # multi-environment main genotypic model
#' Env <- as.factor(c(2,3)) #subset of 2 environments
#' pheno_geno <- data.frame(env = gl(n = 2, k = nrow(Y), labels = Env),
#'                          GID = gl(n = nrow(Y), k = 1,length = nrow(Y) * length(Env)),
#'                          value = as.vector(Y[,2:3]))
#' 
#' K <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MM")
#' y <- pheno_geno[,3]
#' fit <- BGGE(y = y, K = K, ne = rep(nrow(Y), length(Env)), ite = 300, burn = 100,thin = 1)
#' 
#' 
#' @seealso 
#' \code{\link[BGLR]{BGLR}}
#' 
#' @references
#' Cuevas, J., Perez-Elizalde, S., Soberanis, V., Perez-Rodriguez, P., Gianola, D., & Crossa, J. (2014).
#' Bayesian genomic-enabled prediction as an inverse problem. G3: Genes, Genomes, Genetics, 4(10), 1991-2001.
#' 
#' Sousa, M. B., Cuevas, J., Oliveira, E. G. C., Perez-Rodriguez, P., Jarquin, D., Fritsche-Neto, R., Burgueno, J. 
#' & Crossa, J. (2017). Genomic-enabled prediction in maize using kernel models with genotype x environment interaction.
#' G3: Genes, Genomes, Genetics, 7(6), 1995-2014.
#' 
#' 
#' @export
BGGE <- function(y, K, XF = NULL, ne, ite = 1000, burn = 200, thin = 3, verbose = FALSE, tol = 1e-10, R2 = 0.5) {
  
  ### PART I  - Conditional distributions functions and eigen descomposition ####
  # Conditional Distribution of tranformed genetic effects b (U'u)
  #' @import stats
  dcondb <- function(n, media, vari) {
    sd <- sqrt(vari)
    return(rnorm(n, media, sd))
  }
  
  # Conditional distribution of compound variance of genetic effects (sigb)
  dcondsigb <- function(b, deltav, n, nu, Sc) {
    z <- sum(b ^ 2 * deltav)
    return(1 / rgamma(1, (n + nu) / 2, (z + nu * Sc) / 2))
  } 
  
  # Conditional distribution of residual compound variance 
  dcondsigsq <- function(Aux, n, nu, Sce) {
    return(1 / rgamma(1, (n + nu) / 2, crossprod(Aux) / 2 + Sce / 2))
  }
  
  # Conditional fixed effects
  rmvnor <- function(n,media,sigma){
    z <- rnorm(n)
    return( media + crossprod(chol(sigma), z))
  }
  
  # Function for eigen descompositions
  eig <- function(K, tol) {
    ei <- eigen(K)
    fil <- which(ei$values > tol)
    return(list(ei$values[fil], ei$vectors[, fil]))
  }
  
  # Set spectral decomposition
  setDEC <- function(K, tol, ne) {
    
    sK <- vector("list", length = length(K))
    typeM <- sapply(K, function(x) x$Type)
    
    if (!all(typeM %in% c("BD", "D")))
      stop("Matrix should be of types BD or D")
    
    if (missing(ne)) {
      if (any(typeM == "BD"))
        stop("For type BD, number of subjects in each sub matrices should be provided")
    } else {
      if (length(ne) <= 1 & any(typeM == "BD"))
        stop("ne invalid. For type BD, number of subjects in each sub matrices should be provided")
      
      nsubK <- length(ne)
      
      if (nsubK > 1) {
        posf <- cumsum(ne)
        posi <- cumsum(c(1,ne[-length(ne)]))
      }
    }
    
    
    for (i in 1:length(K)) {
      if (K[[i]]$Type == "D") {
        tmp <- list()
        ei <- eig(K[[i]]$Kernel, tol)
        tmp$s <- ei[[1]]
        tmp$U <- ei[[2]]
        tmp$tU <- t(ei[[2]])
        tmp$nr <- length(ei[[1]])
        tmp$type <- "D"
        tmp$pos <- NA
        sK[[i]] <- list(tmp)
      }
      
      if (K[[i]]$Type == "BD") {
        cont <- 0
        tmp <- list()
        for (j in 1:nsubK) {
          Ktemp <- K[[i]]$Kernel[(posi[j]:posf[j]), (posi[j]:posf[j])]
          ei <- eig(Ktemp, tol)
          if (length(ei[[1]]) != 0) {
            cont <- cont + 1
            tmp[[cont]] <- list()
            tmp[[cont]]$s <- ei[[1]]
            tmp[[cont]]$U <- ei[[2]]
            tmp[[cont]]$tU <- t(ei[[2]])
            tmp[[cont]]$nr <- length(ei[[1]])
            tmp[[cont]]$type <- "BD"
            tmp[[cont]]$pos <- c(posi[j], posf[j])
          }
        }
        
        if (length(tmp) > 1) {
          sK[[i]] <- tmp
        } else {
          sK[[i]] <- list(tmp[[1]])
        }
      }
    }
    return(sK)
  }
  
  # verbose part I
  if (as.numeric(verbose) != 0) {
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
  
  if (nNa > 0) { 
    y[whichNa] <- mu
  }
  
  # name of each kernel (important to following procedures)
  if (is.null(names(K))) {
    names(K) <- paste("K", seq(length(K)), sep = "")
  }
  
  # initial values of fixed effects 
  if (!is.null(XF)) {
    Bet <- solve(crossprod(XF), crossprod(XF, y))
    nBet <- length(Bet)
    tXX <- solve(crossprod(XF))
  }
  
  # Eigen descomposition for nk Kernels
  nk <- length(K)
  nr <- numeric(length(K))
  typeM <- sapply(K, function(x) x$Type)
  
  if (!all(typeM %in% c("BD", "D")))
    stop("Matrix should be of types BD or D")
  
  if (length(ne) == 1 & any(typeM == "BD"))
    stop("Type BD should be used only for block diagonal matrices")
  
  Ei <- setDEC(K = K, ne = ne, tol = tol)
  
  # Initial values for Monte Carlo Markov Chains (MCMC)
  
  nCum <- sum(seq(1, ite) %% thin == 0)
  
  chain <- vector("list", length = 3)
  names(chain) <- c("mu", "varE", "K")
  chain$varE <- numeric(nCum)
  chain$mu <- numeric(nCum)
  
  chain$K <- vector("list", length = nk)
  names(chain$K) <- names(K)
  chain$K[seq(nk)] <- list(list(varU = numeric(nCum)))
  
  cpred <- vector('list', length = nk)
  names(cpred) <- names(K)
  cpred[seq(nk)] <- list(U = matrix(NA_integer_, nrow = nCum, ncol = n))
  
  
  nu <- 3
  Sce <- (nu + 2) * (1 - R2) * var(y, na.rm = TRUE)
  Sc <- numeric(length(K))
  for (i in 1:length(K)) {
    Sc[i] <- (nu + 2) * R2 * var(y, na.rm = T)/mean(diag(K[[i]]$Kernel))
  }
  tau <- 0.01
  u <- list()
  for (j in 1:nk) {
    u[[j]] <- rnorm(n, 0, 1 / (2 * n))
  }
  sigsq <- var(y)
  #Sc <- rep(0.01, nk)
  sigb <- rep(0.2, nk)
  
  temp <- y - mu
  
  if (!is.null(XF)) {
    B.mcmc <- matrix(0, nrow = nCum, ncol = nBet)
    temp <- temp - XF %*% Bet
  }

  temp <- temp - Reduce('+', u)
  nSel <- 0
  i <- 1
  
  ### PART III  Fitted model with training data ####
  
  # Iterations of Gibbs sampler
  while (i <= ite) {
    time.init <- proc.time()[3]
    
    # Conditional of mu
    temp <- temp + mu
    mu <- rnorm(1, mean(temp), sqrt(sigsq/n))
    #mu.mcmc[i] <- mu
    temp <- temp - mu
    
    # Conditional of fixed effects
    if (!is.null(XF)) {
      temp <- temp + XF %*% Bet
      vari <- sigsq * tXX
      media <- tXX %*% crossprod(XF, temp)
      Bet <- rmvnor(nBet, media, vari)
      temp <- temp - XF %*% Bet
    }
    
    # Conditionals x Kernel
    for (j in 1:nk) {
      # Sampling genetic effects
      
      if (typeM[j] == "D") {
        temp <- temp + u[[j]]
        d <- crossprod(Ei[[j]][[1]]$U, temp)
        s <- Ei[[j]][[1]]$s
        deltav <- 1 / s
        lambda <- sigb[j]
        vari <- s * lambda / (1 + s * lambda * tau)
        media <- tau * vari * d
        nr <- Ei[[j]][[1]]$nr
        b <- dcondb(nr, media, vari)
        u[[j]] <- crossprod(Ei[[j]][[1]]$tU ,b)
        temp <- temp - u[[j]]
      }
      
      if (typeM[j] == "BD") {
        nsk <- length(Ei[[j]])
        
        if (length(nsk > 1)) {
          temp <- temp + u[[j]]
          d <- NULL
          s <- NULL
          neiv <- numeric(nsk)
          pos <- matrix(NA, ncol = 2, nrow = nsk)
          for (k in 1:nsk) {
            pos[k,] <- Ei[[j]][[k]]$pos
            d <- c(d, crossprod(Ei[[j]][[k]]$U, temp[pos[k, 1]:pos[k, 2]]))
            neiv[k] <- length(Ei[[j]][[k]]$s)
            s <- c(s, Ei[[j]][[k]]$s)
          }
          deltav <- 1/s
          lambda <- sigb[j]
          vari <- s * lambda / (1 + s * lambda * tau)
          media <- tau*vari*d
          nr <- length(s)
          b <- dcondb(nr, media, vari)
          
          posf <- cumsum(neiv)
          posi <- cumsum(c(1, neiv[-length(neiv)]))
          utmp <- numeric(n)
          for (k in 1:nsk) {
            utmp[pos[k, 1]:pos[k, 2] ] <- crossprod(Ei[[j]][[k]]$tU, b[posi[k]:posf[k] ])
          }
          u[[j]] <- utmp
          temp <- temp - u[[j]]
          
        }else{
          temp <- temp + u[[j]]
          pos <- Ei[[j]]$pos
          d <- crossprod(Ei[[j]][[1]]$U, temp[pos[1]:pos[2]])
          s <- Ei[[j]][[1]]$s
          deltav <- 1/s
          lambda <- sigb[j]
          vari <- s * lambda / (1 + s * lambda * tau)
          media <- tau*vari*d
          nr <- Ei[[j]][[1]]$nr
          b <- dcondb(nr, media, vari)
          utmp <- numeric(n)
          utmp[pos[1]:pos[2]] <-  crossprod(Ei[[j]][[1]]$tU, b)
          u[[j]] <- utmp
          temp <- temp - u[[j]]
        }
      }
      
      # Sampling scale hyperparameters and variance of genetic effects
      sigb[j] <- dcondsigb(b, deltav, nr, nu, Sc[j])
    }
    
    # Sampling residual variance 
    res <- temp
    #Sce <- dcondSc(nu, sigsq)
    sigsq <- dcondsigsq(res, n, nu, Sce)
    tau <- 1 / sigsq
    
    # Predicting missing values
    if (nNa > 0) {
      uhat <- Reduce('+', u)
      
      if (!is.null(XF)) {
        aux <- XF[yNA,] %*% Bet
      }else{
        aux <- 0
      }
      
      y[whichNa] <- aux + mu + uhat[whichNa] + rnorm(n = nNa, sd = sqrt(sigsq))
      temp[whichNa] <- y[whichNa] - uhat[whichNa] - aux - mu
    }
    
    # Separating what is for the chain
    if (i %% thin == 0) {
      nSel <- nSel + 1
      chain$varE[nSel] <- sigsq
      chain$mu[nSel] <- mu
      if (!is.null(XF)) {
        B.mcmc[nSel,] <- Bet
      }
      for (j in 1:nk) {
        cpred[[j]][nSel,] <- u[[j]]
        chain$K[[j]]$varU[nSel] <- sigb[j]
      }
    }
    
    # Verbose 
    if (as.numeric(verbose) != 0 & i %% as.numeric(verbose) == 0) {
      time.end <- proc.time()[3]
      cat("Iter: ", i, "time: ", round(time.end - time.init, 3),"\n")
    }
    i <- i + 1
  }

  
  ###### PART IV  Output ######
  #Sampling
  draw <- seq(ite)[seq(ite) %% thin == 0] > burn
  
  mu.est <- mean(chain$mu[draw])
  yHat <- mu.est
  
  if (!is.null(XF)) {
    B <- colMeans(B.mcmc[draw,])
    yHat <- yHat + XF %*% B
  }
  
  u.est <- sapply(cpred, FUN = function(x) colMeans(x[draw, ]) )
  yHat <- yHat + rowSums(u.est)
  
  out <- list()
  out$yHat <- yHat
  out$varE <- mean(chain$varE[draw])
  out$varE.sd <- sd(chain$varE[draw])
  
  out$K <- vector('list', length = nk)
  names(out$K) <- names(cpred)
  for (i in 1:nk) {
    out$K[[i]]$u <- colMeans(cpred[[i]][draw, ])
    out$K[[i]]$u.sd <- apply(cpred[[i]][draw, ], MARGIN = 2, sd)
    out$K[[i]]$varu <- mean(chain$K[[i]]$varU[draw])
    out$K[[i]]$varu.sd <- sd(chain$K[[i]]$varU[draw])
  }
  
  out$chain <- chain
  out$ite <- ite
  out$burn <- burn
  out$thin <- thin
  out$model <- K$model
  out$kernel <- K$kernel
  out$y <- y
  class(out) <- "BGGE"
  return(out)
}
