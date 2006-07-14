##---------------------------------------------------------------------
## The "arsinh" transformation
##
## Note: the constant -log(2*facs[1]) is added to the transformed data
## in order to achieve h_1(y) \approx log(y) for y\to\infty, that is,
## better comparability to the log transformation.
## It has no effect on the generalized log-ratios.
##--------------------------------------------------------------------
vsnh <- function(y, p, strata) {
  if (!is.matrix(y) || !is.numeric(y))
    stop("vsnh: 'y' must be a numeric matrix.\n")
  
  if (!is.array(p) || !is.numeric(p) || any(is.na(p)))
    stop("'p' must be an array with no NAs.\n")

  if(missing(strata)) {
    strata <- rep(as.integer(1), nrow(y))
  } else {
    if(!is.integer(strata) || !is.vector(strata) || 
       length(strata)!=nrow(y) || any(is.na(strata)))
      stop("'strata' must be an integer vector of length nrow(y) with no NAs.")
  }
  nrstrata <- max(strata)
  
  if(nrstrata==1 && length(dim(p))==2)
    dim(p) <- c(1, dim(p))
  
  if(length(dim(p))!=3 || dim(p)[1]!=nrstrata || dim(p)[2]!=ncol(y) || dim(p)[3]!=2)
    stop("'p' has wrong dimensions.")
  if (any(p[,,2]<=0))
    stop("'p' contains invalid values: factors must be non-negative.")

  hy = .Call("vsn_c", y, as.vector(p), strata, as.integer(2), PACKAGE="vsn")

  ## The old, memory-wasting way:
  ## use the recycling rule: p[strata,,1], p[strata,,2], y, and the result of asinh()
  ## are all n*d-matrices. p[strata,1,2] and the result of log() is an n-vector.
  ## hy   = asinh(p[strata,,1] + p[strata,,2] * y) - log(2*p[strata,1,2])

  dimnames(hy) = dimnames(y)
  return(hy)
}

