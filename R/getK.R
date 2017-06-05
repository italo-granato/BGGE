#' compute kernel matrix for GxE models
#'
#' @usage getK(Y, X, XF=NULL, method=c("GK", "G-BLUP"), h=1, model = c("SM", "MM", "MDs", "MDe", "Cov"))
#'
#' @param Y \code{data.frame} Phenotypic data with three columns. The first column is a \code{factor} for assigned environments,
#' the second column is a \code{factor} for assigned individuals and the third column contains the trait of interest.
#' @param X \code{matrix} Marker matrix with individuals in rows and marker in columns
#' @param method Kernel to be used. Methods implemented are the gaussian kernel \code{GK} and the linear kernel \code{G-BLUP}
#' @param h \code{numeric} Bandwidth parameter to create the Gaussian Kernel (GK) matrix. The default for \code{h} is 1, in case don't be provided.
#' Estimation of h can be made using xxxx function.
#' @param model Specifies the genotype by environment model to be fitted. \code{SM} is the single-environment main genotypic effect model,
#' \code{MM} is the multi-environment main genotypic effect model, \code{MDs} is the multi-environment single variance genotype × environment deviation model,
#' \code{MDe} is the multi-environment, environment-specific variance genotype × environment deviation model
#'
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

#' export
getK <- function(Y, X, method = c("GK", "G-BLUP"), h = 1, model = c("SM", "MM", "MDs", "MDe", "Cov"))
{
  #Force to data.frame
  Y <- as.data.frame(Y)
  Y[colnames(Y)[1:2]] <- lapply(Y[colnames(Y)[1:2]], factor)
  
  subj <- levels(Y[,2])
  env <- levels(Y[,1])
  nEnv <- length(env)
  
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
  
  
  return(list(K = out, y = Y[,3]))#) #FALSE
}
