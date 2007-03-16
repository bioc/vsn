## See further below for class definitions.
## First come the validity functions, because they are needed for that

validScalar = function(ob, nm, min=0, max=+Inf) {
  s = slot(ob, nm)
  if((length(s)!=1)||any(is.na(s))||(s<min)||(s>max))
    stop(sprintf("'%s' must be a scalar with values between %g and %g.", nm, min, max))
}
validLogical = function(ob, nm) {
  s = slot(ob, nm)
  if((length(s)!=1)||any(is.na(s)))
    stop(sprintf("'%s' must be a logical of length one and not be NA.", nm))
}

equalOrZero = function(i, j) ((i==j)||(i==0))

validityVsnInput = function(object){
  validScalar(object, "subsample", min=0, max=nrow(object))
  validScalar(object, "cvg.niter", min=1)
  validScalar(object, "cvg.eps",   min=0)
  validScalar(object, "lts.quantile", min=0.5, max=1)
  validLogical(object, "verbose")
  validLogical(object, "ordered")

  if(!equalOrZero(length(object@strata), nrow(object@x)))
    stop("'length(strata)' must must match 'nrow(x)'.")

  if(any(is.na(object@strata)))
    stop("'strata' must not contain NA values.")

  if(!is.numeric(object@pstart)||length(dim(object@pstart))!=3)
    stop("'pstart' must be a 3D array.")
  
  if(!all(dim(object@pstart)==c(nlevels(object@strata), ncol(object@x), 2)))
    stop("Invalid dimensions of 'pstart'.")

  return(TRUE)
}


## strata may be of length 0 (in which case there are no strata).
## refh may be of length 0 (in which case there is no reference).
## hx may be of length 0 (in which case the data is not provided).
## If any these slots is of length>0, then they must agree in size.
validityVsn = function(object){
  if(any(is.na(object@par))||(length(dim(object@par))!=3))
    stop("'par' must be a 3D array and not contain NA values.")

  if(dim(object@par)[3]!=2)
    stop("'dim(par)[3]' must be equal to 2.")
  
  if((length(object@n)!=1)||any(is.na(object@n)))
    stop("'n' must be of length 1 and not NA.")

  if(any(is.na(object@refh)))
    warning("'refh' contains NA.")
  
  if(length(object@refsigma)!=1)
    stop("'refsigma' must be of length 1.")

  if(!equalOrZero(ncol(object@hx), dim(object@par)[2]))
    stop("'ncol(hx)' and 'dim(object@par)[2]' must match.")

  if(!equalOrZero(length(object@strata), object@n))
    stop("'length(strata)' must match 'n'.")

  if(!equalOrZero(length(object@refh), object@n))
    stop("'length(refh)' must match 'n'.")

  if(!equalOrZero(nrow(object@hx), object@n))
    stop("'nrow(hx)' must match 'n'.")

  if(length(object@strata)>0)
    if(nlevels(object@strata)!=dim(object@par)[1])
      stop("'nlevels(strata)' and 'dim(par)[1]' must match.")
  
  return(TRUE)
}

##------------------------------------------------------------
## Class vsn
##------------------------------------------------------------
setClass("vsn",
  representation(
    par    = "array",
    n      = "integer",
    strata = "factor", 
    refh   = "numeric",
    refsigma  = "numeric",
    hx   = "matrix"),
  prototype = list(
    par    = array(0, dim=c(0,0,2)),
    n      = as.integer(0),
    strata = factor(integer(0), levels="all"),
    refh   = numeric(0),
    refsigma  = as.numeric(NA),
    hx   = matrix(0, nrow=0, ncol=0)),
  validity = validityVsn)

##------------------------------------------------------------
## Class vsnInput
##------------------------------------------------------------

setClass("vsnInput",
  representation(
    x  = "matrix",     ## The n*d data matrix
    reference = "vsn", ## A result from a previous fit (for reference normalization)          
    strata = "factor", ## Factor of length n, aligned with rows of x. Special case:
                       ##  It can also be length 0, in case there is only one stratum
    ordered = "logical",  ## Are the levels consecutive in "strata"?               
    lts.quantile = "numeric",
    subsample = "integer",
    verbose   = "logical",
    pstart    = "array",     ## Start parameters (3D array: nrstrata * d * 2)
    cvg.niter = "integer",
    cvg.eps   = "numeric"),
  prototype = list(
    x = matrix(as.numeric(NA), nrow=0, ncol=0),
    reference = new("vsn"),
    strata = factor(integer(0), levels="all"),
    ordered = FALSE,
    lts.quantile = 1,
    subsample = as.integer(0),
    verbose = TRUE,
    pstart = array(as.numeric(NA), dim=c(1,0,2)),
    cvg.niter = as.integer(1),
    cvg.eps  = 0),          
  validity = validityVsnInput)
                  

