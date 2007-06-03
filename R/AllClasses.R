## See further below for class definitions.
## First come the validity functions, because they are needed for that

validScalarNumericSlot = function(ob, nm, min=0, max=+Inf) {
  ## S4 already makes sure s is numeric 
  s = slot(ob, nm)
  if((length(s)!=1)||any(is.na(s))||(s<min)||(s>max))
    return(sprintf("'%s' must be a numeric vector of length 1 with values between %g and %g.", nm, min, max))
  TRUE
}

validScalarDoubleListElt = function(ob, nm, min=0, max=+Inf) {
  s = ob[[nm]]
  if((!is.double(s))||(length(s)!=1)||any(is.na(s))||(s<min)||(s>max))
    return(sprintf("'%s' must be a numeric double precision vector of length 1 with values between %g and %g.", nm, min, max))
  TRUE
}

validScalarIntListElt = function(ob, nm, min=0L, max=.Machine$integer.max) {
  s = ob[[nm]]
  if((!is.integer(s))||(length(s)!=1)||any(is.na(s))||(s<min)||(s>max))
    return(sprintf("'%s' must be an integer vector of length 1 with values between %g and %g.", nm, min, max))
  TRUE
}

validLogical = function(ob, nm) {
  s = slot(ob, nm)
  if((length(s)!=1)||any(is.na(s)))
    return(sprintf("'%s' must be a logical of length one and not be NA.", nm))
  TRUE
}

equalOrZero = function(i, j) ((i==j)||(i==0))

validityVsnInput = function(object){
  r = validScalarNumericSlot(object, "subsample", min=0, max=nrow(object))
  if(!identical(r, TRUE)) return(r)

  r = validScalarNumericSlot(object, "lts.quantile", min=0.5, max=1)
  if(!identical(r, TRUE)) return(r)

  ## could also define class for this, might be a bit less tedious
  r = validLogical(object, "verbose")
  if(!identical(r, TRUE)) return(r)

  r = validLogical(object, "ordered")
  if(!identical(r, TRUE)) return(r)

  if(!equalOrZero(length(object@strata), nrow(object@x)))
    return("'length(strata)' must must match 'nrow(x)'.")

  if(any(is.na(object@strata)))
    return("'strata' must not contain NA values.")

  if(!is.numeric(object@pstart)||length(dim(object@pstart))!=3)
    return("'pstart' must be a 3D array.")
  
  if(!all(dim(object@pstart)==c(nlevels(object@strata), ncol(object@x), 2)))
    return("Invalid dimensions of 'pstart'.")

  if(!is.list(object@optimpar)||(length(object@optimpar)!=7)||
     !identical(names(object@optimpar), c("factr", "pgtol", "lower", "maxit", "trace", "cvg.niter", "cvg.eps")))
    return("'optimpar' must be a list with elements 'factr', 'pgtol', 'lower', 'maxit', 'trace', 'cvg.niter', 'cvg.eps'.")

  r = validScalarDoubleListElt(object@optimpar, "factr")
  if(!identical(r, TRUE)) return(r)
  r = validScalarDoubleListElt(object@optimpar, "pgtol")
  if(!identical(r, TRUE)) return(r)
  r = validScalarDoubleListElt(object@optimpar, "lower")
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "maxit", min=1)
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "trace")
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "cvg.niter", min=1)
  if(!identical(r, TRUE)) return(r)
  r = validScalarDoubleListElt(object@optimpar, "cvg.eps")
  if(!identical(r, TRUE)) return(r)

  return(TRUE)
}


## strata may be of length 0 (in which case there are no strata).
## refh may be of length 0 (in which case there is no reference).
## hx may be of length 0 (in which case the data is not provided).
## If any these slots is of length>0, then they must agree in size.
validityVsn = function(object){
  if(any(is.na(object@par))||(length(dim(object@par))!=3))
    return("'par' must be a 3D array and not contain NA values.")

  if(dim(object@par)[3]!=2)
    return("'dim(par)[3]' must be equal to 2.")
  
  if((length(object@n)!=1)||any(is.na(object@n)))
    return("'n' must be of length 1 and not NA.")

##  if(any(is.na(object@refh)))
##    warning("'refh' contains NA.")
  
  if(length(object@refsigma)!=1)
    return("'refsigma' must be of length 1.")

  if(!equalOrZero(ncol(object@hx), dim(object@par)[2]))
    return("'ncol(hx)' and 'dim(object@par)[2]' must match.")

  if(!equalOrZero(length(object@strata), object@n))
    return("'length(strata)' must match 'n'.")

  if(!equalOrZero(length(object@refh), object@n))
    return("'length(refh)' must match 'n'.")

  if(!equalOrZero(nrow(object@hx), object@n))
    return("'nrow(hx)' must match 'n'.")

  if(length(object@strata)>0)
    if(nlevels(object@strata)!=dim(object@par)[1])
      return("'nlevels(strata)' and 'dim(par)[1]' must match.")
  
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
    optimpar  = "list"),     ## factr, pgtol, lower, maxit, trace, cvg.niter, cvg.eps
  prototype = list(
    x = matrix(as.numeric(NA), nrow=0, ncol=0),
    reference = new("vsn"),
    strata = factor(integer(0), levels="all"),
    ordered = FALSE,
    lts.quantile = 1,
    subsample = 0L,
    verbose = TRUE,
    pstart = array(as.numeric(NA), dim=c(1,0,2)),
    optimpar = list(factr=5e7, pgtol=0, lower=2e-5,
                 maxit=60000L, trace=0L, cvg.niter=1L, cvg.eps=0)),
  validity = validityVsnInput)
                  

