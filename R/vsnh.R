##---------------------------------------------------------------------
## The "arsinh" transformation
## note: the constant -log(2*facs[1]) is added to the transformed data
## in order to achieve h_1(y) \approx log(y) for y\to\infty, that is,
## better comparability to the log transformation.
## It has no effect on the generalized log-ratios.
##---------------------------------------------------------------------
vsnh = function(y, p) {
  if (!is.matrix(y) || !is.numeric(y))
    stop("vsnh: argument y must be a numeric matrix.\n")
  if (!is.vector(p) || !is.numeric(p) || any(is.na(p)))
    stop("vsnh: argument p must be a numeric vector with no NAs.\n")
  if (2*ncol(y) != length(p))
    stop("vsnh: argument p must be a vector of length 2*ncol(y).\n")

  offs = matrix(p[         1  :  ncol(y) ], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  facs = matrix(p[(ncol(y)+1):(2*ncol(y))], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  hy   = asinh(offs + facs * y) - log(2*facs[1])
  dimnames(hy) = dimnames(y)
  return(hy)
}
