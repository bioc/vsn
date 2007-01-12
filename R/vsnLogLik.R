##----------------------------------------------------------------
## given vsnInput object (containg data matrix and possibly strata)
## calculate likelihood for parameters p
##----------------------------------------------------------------

vsnLogLik = function(object, p, refh, refsigma) {
  if(!missing(refh))
    stopifnot(length(refh)==nrow(object), length(refsigma)==1,
              is.numeric(refsigma), !is.na(refsigma))
  else {
    if(!missing(refsigma))
      stop("Arguments 'refh' and 'refsigma' must be either both present or both missing.")
    refh=numeric(0)
    refsigma=0
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
