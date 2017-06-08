#' Genotype x Environment models using linear or gaussian kernel
#'
#' @usage mainGE(Y, X, XF = NULL, W = NULL, kernel = c("GK", "GB"),
#'               h = 1, model = c("SM", "MM", "MDs", "MDe", "Cov"),
#'               nIter = 1000, burnIn = 300, thin = 5, verbose = FALSE, me = 1e-10,...)
#'
#' @param Y \code{data.frame} Phenotypic data with three columns. The first column is a \code{factor} for assigned environments,
#' the second column is a \code{factor} for assigned individuals and the third column contains the trait of interest.
#' @param X \code{matrix} Marker matrix with individuals in rows and marker in columns
#' @param XF \code{matrix} Design matrix (\eqn{n \times p}) for fixed effects
#' @param W \code{matrix} Environmental covariance matrix of dimension (\eqn{n \times n}). If \code{W} is provided, it can be used along with \code{MM} and \code{MDs} models. See details
#' @param kernel Kernel to be used. Methods implemented are the gaussian kernel \code{Gk-EB} and the linear kernel \code{G-BLUP}
#' @param h \code{numeric} Bandwidth parameter to create the gaussian kernel matrix. For \code{method} \code{Gk-EB}, if \code{h} is not provided,
#' then it is computed following a empirical bayesian selection method. See details
#' @param model Specifies the genotype by environment model to be fitted. \code{SM} is the single-environment main genotypic effect model,
#' \code{MM} is the multi-environment main genotypic effect model, \code{MDs} is the multi-environment single variance genotype × environment deviation model,
#' \code{MDe} is the multi-environment, environment-specific variance genotype × environment deviation model
#' @param nIter \code{integer} Number of iterations.
#' @param burnIn \code{integer} Number of iterations to be discarded as burn-in.
#' @param thin \code{integer} Thinin interval.
#' @param verbose \code{logical} Should report be printed on screen?
#' @param me \code{numeric} tolerance for zero. Default is 1e-10
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
#' \end{itemize}
#'
#' @return
#' return variance components and 
#'
#' @seealso \code{\link[MTM]{MTM}}, \code{\link{getK}} and \code{\link{BLMMD}}
#'
#'@examples
#' # create kernel matrix for model MDs using wheat dataset
#' library(BGLR)
#' 
#' data(wheat)
#' X <- scale(wheat.X, scale = TRUE, center = TRUE)
#' rownames(X) <- 1:599
#' pheno_geno <- data.frame(env = gl(n = 4, k = 599), 
#'                GID = gl(n=599, k=1, length = 599*4),
#'                value = as.vector(wheat.Y))
#' 
#' fit <- mainGE(Y = pheno_geno, X = X, method = GB, model = "MDs")
#'  
#'
#'export
mainGE <- function(Y, X, XF=NULL, W=NULL, kernel=c("GK", "GB"), h=1, model = c("SM", "MM", "MDs", "MDe", "Cov"),
                   nIter = 1000, burnIn = 300, thin = 5, verbose = FALSE, me = 1e-10, ...) {
  hasW <- !is.null(W)
  
  names(Y) <- c("environ", "subjects", "value")
  subj <- levels(Y$subjects)
  env <- levels(Y$environ)
  
  #setting kernels
  setK <- getK(Y = Y, X = X, kernel = kernel, h = h, model = model)
  K <- setK$K
  y <- as.vector(setK$y)
  
  if (model %in% c("MM", "MDs")) {
    if (hasW) {
      Ze <- model.matrix(~factor(Y[,1])-1)
      Zg <- model.matrix(~factor(Y[,2])-1)
      EC <- Ze %*% (W %*% t(Ze))
      GEC <- setK$K$G * EC
      K <- c(K, list(EC = list(K = EC), GEC = list(K = GEC)))
    }
  }
  
  if (model == "Cov") {
    tmpY <- matrix(nrow = length(subj), ncol = length(env), data = NA)
    colnames(tmpY) <- env
    rownames(tmpY) <- subj
    
    for (i in 1:length(env)) {
      curEnv <- Y$env == env[i]
      curSub <- match(Y[curEnv, 2], rownames(tmpY))
      tmpY[curSub, i] <- Y[curEnv, 3]
    }
    
    fit <- GEcov(Y = tmpY, K = K, nIter = nIter, burnIn = burnIn, thin = thin,...)
  }
  else {
     fit <- BLMMD(y = y, K = K, XF = XF, ite = nIter, burn = burnIn, thin = thin, verbose = verbose, me=me)
  }
  return(fit)
}


GEcov <- function(Y, K, nIter = 1000, burnIn = 300, thin = 5, ...)
{
  env <- ncol(Y)
  In <- diag(x = 1, nrow = nrow(Y), ncol = nrow(Y))
  
  # How to find df0 and S0?
  ETA <- list(list( K = K, COV = list(type = 'UN', df0 =  env, S0 = diag(env))),
              list(K = In, COV = list(type = 'UN', df0 = env, S0 = diag(env))))
  resCov <- list( type = 'DIAG', S0 = rep(1, env), df0 = rep(1, env))
  
  #' @importFrom MTM MTM
  fit <- MTM(Y = Y, K = ETA, resCov = resCov, nIter = nIter, burnIn = burnIn, thin = thin,...)
  
  return(fit)
}