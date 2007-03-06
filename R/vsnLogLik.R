##----------------------------------------------------------------
## given vsnInput object (containg data matrix and possibly strata)
## calculate likelihood for parameters p
##----------------------------------------------------------------

vsnLogLik = function(object, p, refh=numeric(0), refsigma=as.numeric(NA)) {
  if(length(refh)>0) {
    if(length(refh)!=nrow(object))
      stop("length(refh) must be either 0 or nrow(object).")
    if((length(refsigma)!=1) || !is.numeric(refsigma) || !is.finite(refsigma))
      stop("'refsigma' must be a finite scalar numeric not NA.")
  } else {
    if(!is.na(refsigma))
      stop("If length(refh) is 0, refsigma must not be specified.")
  }

  if(is.vector(p))
    dim(p)=c(length(p), 1)  
    
  res = matrix(as.numeric(NA), nrow=1+nrow(p), ncol=ncol(p))
  istrat = calcistrat(object) 
  for(j in 1:ncol(p))
    res[, j] = .Call("vsn2_point", object@x, p[,j], istrat, refh, refsigma, PACKAGE="vsn")

  ## invert sign since the above calculate the _negative_ log likelihood
  return(-res)
}
