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
