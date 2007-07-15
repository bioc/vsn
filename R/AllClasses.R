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
    return(sprintf("'%s' must be an integer vector of length 1 with values between %d and %d.", nm, min, max))
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

  compulsoryElements = c("factr", "pgtol", "lower", "maxit", "trace", "REPORT", "cvg.niter", "cvg.eps")
  if(!is.list(object@optimpar)||!identical(names(object@optimpar), compulsoryElements))
    return(paste("'optimpar' must be a list with elements ",
                 paste("'", compulsoryElements, "'", collapse=", ", sep=""), ".", sep=""))

  r = validScalarDoubleListElt(object@optimpar, "factr")
  if(!identical(r, TRUE)) return(r)
  r = validScalarDoubleListElt(object@optimpar, "pgtol")
  if(!identical(r, TRUE)) return(r)
  r = validScalarDoubleListElt(object@optimpar, "lower")
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "maxit", min=1L)
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "trace", min=0L, max=6L)
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "REPORT")
  if(!identical(r, TRUE)) return(r)
  r = validScalarIntListElt(object@optimpar, "cvg.niter", min=1L)
  if(!identical(r, TRUE)) return(r)
  r = validScalarDoubleListElt(object@optimpar, "cvg.eps")
  if(!identical(r, TRUE)) return(r)

  return(TRUE)
}


## strata may be of length 0 (in which case there are no strata).
## hx may have 0 rows (in which case the transformed data has not (yet) been computed).
## If length(strata) or nrow(hx) are >0, then they must be the same
##   as length(mu) (which is always >0).
validityVsn = function(object){
  if(any(is.na(object@coefficients))||(length(dim(object@coefficients))!=3))
    return("'coefficients' must be a 3D array and not contain NA values.")

  if(dim(object@coefficients)[3]!=2)
    return("'dim(coefficients)[3]' must be equal to 2.")
  
  if(length(object@sigsq)!=1)
    return("'sigsq' must be of length 1.")

  if(!equalOrZero(ncol(object@hx), dim(object@coefficients)[2]))
    return("'ncol(hx)' and 'dim(object@coefficients)[2]' must match.")

  if(!equalOrZero(length(object@strata), length(object@mu)))
    return("'length(strata)' must be 0 or equal to 'length(mu)'.")

  if(!equalOrZero(nrow(object@hx), length(object@mu)))
    return("'nrow(hx)' must be 0 or equal to 'length(mu)'.")

  if(length(object@strata)>0)
    if(nlevels(object@strata)!=dim(object@coefficients)[1])
      return("'nlevels(strata)' and 'dim(coefficients)[1]' must match.")

  if(!((length(object@lbfgsb)==1)&&(is.integer(object@lbfgsb))))
    return("'lbfgsb' must be an integer of length 1.")
  
  return(TRUE)
}

##------------------------------------------------------------------------------------
## Class vsn
## Slot 'coefficients' is a 3D array: nrstrata * d * 2, with 2 parameters for each stratum and
##   array. The first of these 2 is the background offset off, the second the scaling
##   factor fac. The transformation is hence y -> asinh((y+off)/fac)
##------------------------------------------------------------------------------------
setClass("vsn",
  representation(
    coefficients = "array",
    strata = "factor", 
    mu     = "numeric",
    sigsq  = "numeric",
    hx     = "matrix",
    lbfgsb = "integer"),
  prototype = list(
    coefficients = array(0, dim=c(0,0,2)),
    strata       = factor(integer(0), levels="all"),
    mu           = numeric(0),
    sigsq        = as.numeric(NA),
    hx           = matrix(0, nrow=0, ncol=0),
    lbfgsb       = as.integer(NA)),
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
    pstart    = "array",     ## Start parameters: see comment on slot 'coefficients' in definition of class 'vsn'
    optimpar  = "list"),     
  prototype = list(
    x = matrix(as.numeric(NA), nrow=0, ncol=0),
    reference = new("vsn"),
    strata = factor(integer(0), levels="all"),
    ordered = FALSE,
    lts.quantile = 1,
    subsample = 0L,
    verbose = TRUE,
    pstart = array(as.numeric(NA), dim=c(1L,0L,2L)),
    optimpar = list(factr=5e7, pgtol=2e-4, lower=2e-4,
                 maxit=60000L, trace=0L, REPORT=10L, cvg.niter=7L, cvg.eps=0)),
  validity = validityVsnInput)
                  

