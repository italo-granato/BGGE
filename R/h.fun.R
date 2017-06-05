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