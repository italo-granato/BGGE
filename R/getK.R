#' Kernel matrix for GE genomic selection models
#' 
#' Create kernel matrix for GE genomic prediction models 
#'
#' @usage getK(Y, X, kernel = c("GK", "GB"), K = NULL, h = 1,
#'              model = c("SM", "MM", "MDs", "MDe"), quantil = 0.5,
#'              geno.as.random = FALSE)
#'
#' @param Y \code{data.frame} Phenotypic data with three columns. The first column is a factor for environments,
#' the second column is a factor identifying genotypes, and the third column contains the trait of interest
#' @param X Marker matrix with individuals in rows and markers in columns. Missing markers are not allowed.
#' @param kernel Kernel to be created internally. Methods currently implemented are the Gaussian \code{GK} and the linear \code{GBLUP} kernel
#' @param K \code{matrix} Single kernel matrix in case it is necessary to use a different kernel from \code{GK} or \code{GBLUP}
#' @param h \code{numeric} Bandwidth parameter to create the Gaussian Kernel (GK) matrix. The default for \code{h} is 1.
#' Estimation of h can be made using a Bayesian approach as presented in Perez-Elizalde et al. (2015)
#' @param model Specifies the genotype \eqn{\times} environment model to be fitted. It currently supported the 
#' models  \code{SM}, \code{MM}, \code{MDs} and \code{MDe}. See Details
#' @param quantil Specifies the quantile to create the Gaussian kernel.
#' @param geno.as.random if \code{TRUE}, kernel related to random intercept of genotype is included.
#' 
#' @details
#' The aim is to create kernels to fit GE interaction models applied to genomic prediction.
#' Two standard genomic kernels are currently supported:
#' \code{GB} creates a linear kernel resulted from the cross-product of centered and standardized 
#' marker genotypes divide by the number of markers \eqn{p}:
#'     \deqn{GB = \frac{XX^T}{p}}
#' Another alternative is the Gaussian Kernel \code{GK}, resulted from:
#'  \deqn{ GK (x_i, x_{i'}) = exp(\frac{-h d_{ii'}^2}{q(d)})}
#' where \eqn{d_{ii'}^2} is the genetic distance between individuals based on markers scaled 
#' by some percentile \eqn{{q(d)}} and \eqn{h} is the bandwidth parameter. However, 
#' other kernels can be provided through \code{K}. In this case, arguments \code{X}, 
#' \code{kernel} and \code{h} are ignored.
#' 
#' Currently, the supported models for GE kernels are:
#' \itemize{
#' \item \code{SM}: is the single-environment main genotypic effect model - It fits the data for a 
#' single environment, and only one kernel is produced. 
#' \item \code{MM}: is the multi-environment main genotypic effect model - It consideres the main
#' random genetic effects across environments. Thus, just one kernel is produced, of order 
#' \eqn{n \times n}, related to the main effect across environments. 
#' \item \code{MDs}: is the multi-environment single variance genotype x environment deviation 
#' model - It is an extension of \code{MM} by adding the random interaction effect of 
#' environments with genotype information. Thus, two kernels are created, one related to the 
#' main effect across environment, and the second is associated with single genotype by environment effect.
#' \item \code{MDe}: is the multi-environment, environment-specific variance genotype x environment 
#' deviation model - It separates the genetic effects into the main genetic 
#' effects and the specific genetic effects (for each environment). Thus, one kernel 
#' for across environments effect and \eqn{j} kernels are created, one for each 
#' environment.
#' }
#' These GE genomic models were compared and named by Sousa et al. (2017) and can be increased by using 
#' the kernel related to random intercept of genotype through \code{geno.as.random}.
#' 
#' @return
#' This function returns a two-level list, which specifies the kernel and the type of matrix. 
#' The latter is a classification according to its structure, i. e.,
#'  if the matrix is dense or a block diagonal. For the main effect (\code{G}), 
#'  the matrix is classified as dense (D). On the other hand, matrices for environment-specific and 
#'  genotype by environment effect (\code{GE}) are considered diagonal block (BD). This classification is used 
#'  as part of the prediction through the BGGE function.
#'  
#' @references
#' Jarquin, D., J. Crossa, X. Lacaze, P. Du Cheyron, J. Daucourt, J. Lorgeou, F. Piraux, L. Guerreiro, P. Pérez, M. Calus, J. Burgueño,
#' and G. de los Campos. 2014. A reaction norm model for genomic selection using high-dimensional genomic and 
#' environmental data. Theor. Appl. Genet. 127(3): 595-607.
#' 
#' Lopez-Cruz, M., J. Crossa, D. Bonnett, S. Dreisigacker, J. Poland, J.-L. Jannink, R.P. Singh, E. Autrique,
#' and G. de los Campos. 2015. Increased prediction accuracy in wheat breeding trials using a marker × environment
#' interaction genomic selection model. G3: Genes, Genomes, Genetics. 5(4): 569-82.
#' 
#' Perez- Elizalde, S. J. Cuevas, P. Perez-Rodriguez, and J. Crossa. 2015. Selection of the
#' Bandwidth Parameter in a Bayesian Kernel Regression Model for Genomic-Enabled Prediction. 
#' Journal of Agricultural, Biological, and Environmental Statistics (JABES), 20(4):512-532.
#' 
#' Sousa, M. B., Cuevas, J., Oliveira, E. G. C., Perez-Rodriguez, P., Jarquin, D., Fritsche-Neto, R., Burgueno, J. 
#' & Crossa, J. (2017). Genomic-enabled prediction in maize using kernel models with genotype x environment interaction.
#' G3: Genes, Genomes, Genetics, 7(6), 1995-2014.
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
#'                
#'  K <- getK(Y = pheno_geno, X = X, kernel = "GB", model = "MDs")              
#' 
#' 
#' 
#' @export
getK <- function(Y, X, kernel = c("GK", "GB"), K = NULL, h = 1, model = c("SM", "MM", "MDs", "MDe"), quantil = 0.5,
                 geno.as.random = FALSE)
{
  #Force to data.frame
  Y <- as.data.frame(Y)
  Y[colnames(Y)[1:2]] <- lapply(Y[colnames(Y)[1:2]], factor)
  
  subjects <- levels(Y[,2])
  env <- levels(Y[,1])
  nEnv <- length(env)
  
  switch(model,
         'SM' = {
           if (nEnv > 1)
             stop("Single model choosen, but more than one environment is in the phenotype file")
           Zg <- model.matrix(~factor(Y[,2L]) - 1)
         },
         'Cov' = {
           Zg <- model.matrix(~factor(subjects) - 1)
         },{
           Ze <- model.matrix(~factor(Y[,1L]) - 1)
           Zg <- model.matrix(~factor(Y[,2L]) - 1)
         })
  
  if(is.null(K)){
    if(is.null(rownames(X)))
      stop("Genotype names are missing")
    
    if(!all(subjects %in% rownames(X)))
      stop("Not all genotypes presents in the phenotypic file are in marker matrix")
    
    X <- X[subjects,]
    
    switch(kernel,
           'GB' = {
             # case 'G-BLUP'...
             ker.tmp <- tcrossprod(X) / ncol(X)
             #G <- list(Zg %*% tcrossprod(ker.tmp, Zg))
             G <- list(list(Kernel = Zg %*% tcrossprod(ker.tmp, Zg), Type = "D"))
           },
           'GK' = {
             # case 'GK'...
             D <- (as.matrix(dist(X))) ^ 2
             
             G <- list()
             for(i in 1:length(h)){
               ker.tmp <- exp(-h[i] * D / quantile(D, quantil))
               #G[[i]] <- Zg %*% tcrossprod(ker.tmp, Zg)
               G[[i]] <- list(Kernel = Zg %*% tcrossprod(ker.tmp, Zg), Type = "D")
               }
             },
           {
             stop("kernel selected is not available. Please choose one method available or make available other kernel through argument K")
           })
    }else{
      if(is.null(rownames(K)) | is.null(colnames(K)))
        stop("Genotype names are missing")
      
      # Condition to check if all genotypes name are compatible
      if(!all(subjects %in% rownames(K)))
        stop("Not all genotypes presents in phenotypic file are in the kernel matrix")
      
      K <- K[subjects, subjects] # reordering kernel
      
      ker.tmp <- K
      #G <- list(Zg %*% tcrossprod(ker.tmp, Zg))
      G <- list(list(Kernel = Zg %*% tcrossprod(ker.tmp, Zg), Type = "D"))
  }
  
    names(G) <- if(length(G) > 1) paste("G", seq(length(G)), sep ="_") else "G"

  switch(model,
       
         'SM' = {
           out <- G
         }, 
         
         'MM' = {
           out <- G
         },
         
         'MDs' = {
           E <- tcrossprod(Ze)
           #GE <- Map('*', G, list(E))
           GE <- lapply(G, function(x) list(Kernel = x$Kernel * E, Type = "BD"))
           names(GE) <- if(length(G) > 1) paste("GE", seq(length(G)), sep ="_") else "GE"
           out <- c(G, GE)
         },
         
         'MDe' = {
           ZEE <- matrix(data = 0, nrow = nrow(Ze), ncol = ncol(Ze))
           
           out.tmp <- list()
           
           for(j in 1:length(G)){
           out.tmp <- c(out.tmp, lapply(1:nEnv, function(i){
             ZEE[,i] <- Ze[,i]
             ZEEZ <- ZEE %*% t(Ze)
             #K3 <- G[[j]] * ZEEZ
             K3 <- list(Kernel = G[[j]]$Kernel * ZEEZ, Type = "BD")
             return(K3)
           }))
           }
           name.tmp <- paste(rep(env, length(G)), rep(1:length(G), each = nEnv), sep = "_")
           names(out.tmp) <- if(length(G) > 1) name.tmp else env
           out <- c(G, out.tmp)
           }, #DEFAULT CASE
         {
           stop("Model selected is not available ")
         })
    
    if(geno.as.random){
      Gi <- list(Kernel = Zg %*% tcrossprod(diag(length(subjects)), Zg), Type = "D")
      out <- c(out, list(Gi = Gi))
    }
  
  return(out)
}
