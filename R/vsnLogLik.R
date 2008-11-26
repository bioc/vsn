##----------------------------------------------------------------
## given vsnInput object (containg data matrix and possibly strata)
## calculate likelihood for parameters p
##----------------------------------------------------------------
checkArgs = function(object, mu, sigsq) {

  if(length(mu)>0) {
    if(length(mu)!=nrow(object))
      stop("length(mu) must be either 0 or nrow(object).")
    if((length(sigsq)!=1) || !is.numeric(sigsq) || !is.finite(sigsq))
      stop("'sigsq' must be a finite scalar numeric not NA.")
  } else {
    if(!is.na(sigsq))
      stop("If length(mu) is 0, sigsq must not be specified.")
  }
}

vsnLogLik = function(object, p, mu=numeric(0), sigsq=as.numeric(NA), calib="affine") {
  
  checkArgs(object, mu, sigsq)
  if(is.vector(p))
    dim(p)=c(length(p), 1)  
    
  res = matrix(as.numeric(NA), nrow=1+nrow(p), ncol=ncol(p))
  istrat = calcistrat(object) 
  for(j in 1:ncol(p))
    res[, j] = .Call("vsn2_point",
         object@x,
         p[,j],
         istrat,
         mu,
         sigsq,
         calibCharToInt(calib),
         PACKAGE="vsn")

  ## invert sign since the above calculates the _negative_ log likelihood
  return(-res)
}

vsnHessian = function(object, p, mu=numeric(0), sigsq=as.numeric(NA), eps=1e-4) {
  checkArgs(object, mu, sigsq)

  if((length(dim(p))!=3)||(any(dim(p)!=c(nlevels(object@strata), ncol(object), 2L))))
    stop("'p' has wrong dimensions.")

  istrat = calcistrat(object)

  np = length(p)
  pp = array(p, dim=c(np, np, 2L))
  dp = ifelse(1:np <= prod(dim(p)[1:2]), eps, p*exp(eps))
  for(j in seq_len(np)) { 
    pp[j, j, 1] = pp[j, j, 1] - dp[j]/2
    pp[j, j, 2] = pp[j, j, 2] + dp[j]/2
  }
  dim(pp) = c(np, np*2L)
  ll = vsnLogLik(object, pp, mu, sigsq)[-1, ]
  dim(ll) = c(np, np, 2L)

  res = matrix(as.numeric(NA), nrow=np, ncol=np)
  for(i in seq_len(np))
    for(j in seq_len(np))
      res[i, j] = (ll[i, j, 2] - ll[i, j, 1])/dp[j]
  
  return(res)
}
