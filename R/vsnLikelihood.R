## given vsnInput object (containg data matrix and possibly reference)
##  calculate likelihood for parameters p

vsnLikelihood = function(v, p, refh, refsigma) {
  if(!missing(refh))
    stopifnot(length(refh)==nrow(v), length(refsigma)==1,
              is.numeric(refsigma), !is.na(refsigma))
  else {
    if(!missing(refsigma))
      stop("Arguments 'refh' and 'refsigma' must be either both present or both missing.")
    refh=numeric(0)
    refsigma=0
  }

  if(is.vector(p)) dim(p)=c(length(p), 1)  
  
  res = matrix(as.numeric(NA), nrow=1+nrow(p), ncol=ncol(p))
  istrat = calcistrat(v) 
  for(j in 1:ncol(p))
    res[, j] = .Call("vsn2_point", v@x, p[,j], istrat, refh, refsigma, PACKAGE="vsn")
  return(res)
}
