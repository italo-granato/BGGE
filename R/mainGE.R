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
#'
mainGE <- function(Y, X, XF=NULL, W=NULL, method=c("GK", "G-BLUP"), h=NULL, model = c("SM", "MM", "MDs", "MDe", "Cov"),
                   nIter = 1000, burnIn = 300, thin = 5, ...) {
  hasW <- !is.null(W)
  subjects <- levels(Y[,2])
  environ <- levels(Y[,1])

  ETA.tmp <- getK(Y=Y, X=X, method = method, h=h, model = model)

  if (model %in% c("MM", "MDs")) {
    if (hasW) {
      Ze <- model.matrix(~factor(Y[,1])-1)
      Zg <- model.matrix(~factor(Y[,2])-1)
      EC <- Ze %*% (W %*% t(Ze))
      GEC <- ETA$G * ENVc
      ETA.tmp <- c(ETA.tmp, list(EC = list(K = EC, model = "RKHS"), GEC = list(K = GEC, model = "RKHS")))
    }
  }

  if (model == "Cov") {
    tmpY <- matrix(nrow = length(subjects), ncol = length(environ), data = NA)
    colnames(tmpY) <- subjects
    rownames(tmpY) <- environ

    for (i in 1:length(environ)) {
      curEnv <- Y[, 1L] == environ[i]
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



getK <- function(Y, X, XF = NULL, method = c("GK", "G-BLUP"), h = NULL, model = c("SM", "MM", "MDs", "MDe"))
getK <- function(Y, X, XF = NULL, method = c("GK", "G-BLUP"), h = NULL, model = c("SM", "MM", "MDs", "MDe", "Cov"))
  {

  hasXF <- !is.null(XF)
  hash <- !is.null(h)
  nEnv <- nlevels(Y[,1L])
  Env <- levels(Y[,1L])
  subj <- levels(Y[,2L])

  if(is.null(rownames(X)))
    stop("Marker name is missing")

  if(!all(levels(Y[,2]) %in% rownames(X)))
    stop("Not all genotypes presents in phenotypic file are in marker matrix")

  X <- X[rownames(X) %in% subj,]

  Ze <- model.matrix( ~ factor(Y[,1L]) - 1)
  Zg <- model.matrix( ~ factor(Y[,2L]) - 1)

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

           if (!hash) {
             source(paste(getwd(), "/margh.fun.R", sep = "")) # calling margh.fun
             sol <- optim(c(1, 1), margh.fun, y = Y0, D = D00, q = quantile(D00, 0.05),
                          method = "L-BFGS-B", lower = c(0.05, 0.05), upper = c(6, 30))
             h <- sol$par[1]
           }

           ker.tmp <- exp(-h * D / quantile(D, 0.05))
           G <- Zg %*% (ker.tmp %*% t(Zg))
         },
         {
           stop("Error, the method selected is not available ")
         })

  switch(model,

    SM = {
      if (nEnv > 1)
        stop("single model choosen, but more than one environment is in the phenotype file")
      out <- list(K = ker.tmp, model = "RKHS")
    },

    MM = {
      out <- list(G = list(K = G, model = "RKHS"))
    },

    MDs = {
      E <- tcrossprod(Ze)
      GE <- G * E
      out <- list(G = list(K = G, model = "RKHS"), GE = list(K = GE, model = "RKHS"))
    },

    MDe = {
      if (method == "G-BLUP") {
        ZEE <- matrix(data = 0, nrow = nrow(Ze),ncol = ncol(Ze))

        out <- lapply(1:nEnv, function(i){
          ZEE[,i] <- Ze[,i]
          ZEEZ <- ZEE %*% t(Ze)
          K3 <- G * ZEEZ
          return(list(K = K3, model = "RKHS"))
        })
        out <- c(list(list(K = G, model = "RKHS")), out)
        names(out) <- c("mu", paste("e", 2:length(out), sep = ""))
      }

      else {
        D1 <- list(D00)

        for (i in 1:nEnv) {
          D1[[i + 1L]] <- D[rownames(D) %in% Y[Y[,1L] == Env[i], 2L], rownames(D) %in% Y[Y[,1L] == Env[i], 2L]]

          if (!hash) {
            Y1 <- Y[naY & Y[,1L] == Env[i], ]
            D2 <- D[match(Y1[, 2L], rownames(D)), match(Y1[, 2L], rownames(D))]
            q05 <- quantile(D2, 0.05)
            sol <- optim(c(1, 1), margh.fun, y = Y1[,3L], D = D2, q = q05,
                        method = "L-BFGS-B", lower = c(0.05, 0.05), upper = c(6, 30))
            h[i + 1] <- sol$par[1]
          }
        }

        GKe <- lapply(X = Map("*", D1, -h), FUN = function(x) exp(x/quantile(x, 0.05)))
        v <- integer(nEnv)
        mDiag <- lapply(X = 1:nEnv, FUN = function(i){v[i] <- 1; diag(v)})

        ETA.tmp <- Map(kronecker, mDiag, GKe[2:length(GKe)])
        out <- lapply(c(list(GKe[[1L]]), ETA.tmp), function(x) list(K = x, model = "RKHS"))
      }
    }, #DEFAULT CASE
    {
      stop("Error, the model selected is not available ")
    })

  ifelse(hasXF, #Test
         return(c(list(list(X = XF, model = "FIXED")), out)), #TRUE
         return(out)) #FALSE
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

