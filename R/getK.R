#' compute kernel matrix for GxE models
#'
#' @usage getK(Y, X, kernel = c("GK", "GB"), K = NULL,
#'             h = 1, model = c("SM", "MM", "MDs", "MDe", "Cov"))
#'
#' @param Y \code{data.frame} Phenotypic data with three columns. The first column is a factor for environments,
#' the second column is a factor for genotype identification and the third column contains the trait of interest
#' @param X Marker matrix with individuals in rows and marker in columns
#' @param kernel Kernel to be used. Methods implemented are the gaussian kernel \code{GK} and the linear kernel \code{G-BLUP}
#' @param K Single kernel matrix in case it needs to use a different kernel from those supported
#' @param h \code{numeric} Bandwidth parameter to create the Gaussian Kernel (GK) matrix. The default for \code{h} is 1.
#' Estimation of h can be made using xxxx function
#' @param model Specifies the genotype by environment model to be fitted. It currently supported the models  \code{SM}, \code{MM}, \code{MDs} and \code{MDe}. See Details
#'
#' @details
#' The goal is to create kernels used to fit GxE interaction models. These models can be adjusted using different kernels.
#' \code{GB} creates a linear kernel resulted from the cross-product of centered and standarized marker genotypes divide by the number of markers \eqn{p}:
#'     \deqn{G = \frac{XX^T}{p}}
#' Other kernel option currently supported is the Gaussian Kernel \code{Gk}, resulted from exponential of a genetic distance matrix based on markers scaled by its fifth percentile multiplied by the bandwidth parameter \eqn{h}
#' which has an objective of controlling the rate of decay for correlation between individuals. Thus:
#'  \deqn{ K (x_i, x_{i'}) = exp(-h d_{ii'}^2)}
#' However other kernels can be provided through \code{K}. In this case, arguments \code{X}, \code{kernel} and \code{h} are ignored.
#' 
#' The currently supported models for GxE kernels are:
#' \itemize{
#' \item \code{SM}: is the single-environment main genotypic effect model - The SM model fits the data for single environment
#' \item \code{MM}: is the multi-environment main genotypic effect model - Multi-environment model considering the main random genetic effects across environments.
#' \item \code{MDs}: is the multi-environment single variance genotype x environment deviation model - This model is an extension of \code{MM} by adding the random interaction effect of the environments
#' with the genetic information of genotypes.
#' \item \code{MDe}: is the multi-environment, environment-specific variance genotype x environment deviation model - This model separates the genetic effects of markers into the main marker effects (across environments) and the specific marker effects (for each environment).
#' \item \code{Cov}: is the multi-environment, variance-covariance environment-specific model -
#' }
#' 
#' @return
#' It returns a two-level list, which the first one has the kernels for respective model and the second element is the phenotypic value.
#' 
#' @examples 
#' # create kernel matrix for model MDs using wheat dataset
#' library(BGLR)
#' 
#' data(wheat)
#' X <- scale(wheat.X, scale = TRUE, center = TRUE)
#' rownames(X) <- 1:599
#' pheno_geno <- data.frame(env = gl(n = 4, k = 599), 
#'                GID = gl(n=599, k=1, length = 599*4),
#'                value = as.vector(wheat.Y))
#'  K <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MDs")              
#' 
#' 
#' 
#' @export
getK <- function(Y, X, kernel = c("GK", "GB"), K=NULL, h = 1, model = c("SM", "MM", "MDs", "MDe", "Cov"))
{
  #Force to data.frame
  Y <- as.data.frame(Y)
  Y[colnames(Y)[1:2]] <- lapply(Y[colnames(Y)[1:2]], factor)
  
  subj <- levels(Y[,2])
  env <- levels(Y[,1])
  nEnv <- length(env)
  
  
  if(model == "SM"){
    if (nEnv > 1)
      stop("Single model choosen, but more than one environment is in the phenotype file")
    
    Zg <- model.matrix( ~ factor(Y[,2L]) - 1)
  }else{
    Ze <- model.matrix( ~ factor(Y[,1L]) - 1)
    Zg <- model.matrix( ~ factor(Y[,2L]) - 1)
  }

  if (model == "Cov")
    Zg <- model.matrix( ~ factor(subj) - 1)
  
  if(is.null(K)){
    if(is.null(rownames(X)))
      stop("Genotype names are missing")
    
    if(!all(subj %in% rownames(X)))
      stop("Not all genotypes presents in phenotypic file are in marker matrix")
    
    X <- X[subj,]
    
    switch(kernel,
           'GB' = {
             # case 'G-BLUP'...
             ker.tmp <- tcrossprod(X) / ncol(X)
             G <- Zg %*% tcrossprod(ker.tmp, Zg)
           },
           GK = {
             print("Gk")
             # case 'GK'...
             D <- (as.matrix(dist(X))) ^ 2
             ker.tmp <- exp(-h * D / quantile(D, 0.05))
             G <- Zg %*% tcrossprod(ker.tmp, Zg)
           },
           {
             stop("Method selected is not available ")
           })
    }else{
      # Condition to check if is symmetrical
      if(!all(subj %in% rownames(K)))
        stop("Not all genotypes presents in phenotypic file are in the kernel matrix")
      
      ker.tmp <- K
      G <- Zg %*% tcrossprod(ker.tmp, Zg)
  }
  
  switch(model,
         
         'SM' = {
           out <- list(K = ker.tmp)
         },
         
         'MM' = {
           out <- list(G = list(K = G))
         },
         
         'MDs' = {
           E <- tcrossprod(Ze)
           GE <- G * E
           out <- list(G = list(K = G), GE = list(K = GE))
         },
         
         'MDe' = {
           ZEE <- matrix(data = 0, nrow = nrow(Ze), ncol = ncol(Ze))
           
           out.tmp <- lapply(1:nEnv, function(i){
             ZEE[,i] <- Ze[,i]
             ZEEZ <- ZEE %*% t(Ze)
             K3 <- G * ZEEZ
             return(list(K = K3))
           })
           out <- c(list(list(K = G)), out.tmp)
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

