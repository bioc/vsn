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

compulsoryElements = c("factr", "pgtol", "maxit", "trace", "cvg.niter", "cvg.eps")

calibCharToInt = function(x)
  switch(EXPR=x, affine=0L, none=1L, stop(sprintf("Invalid value '%s' of 'calib'.", x)))


##--------------------------------------------------
## validity method for 'vsnInput'
##--------------------------------------------------
validityVsnInput = function(object){
  if(any(is.nan(object@x)))
    return("The values in 'x' must be finite numeric or NA; please no NaN.")
  
  r = validScalarNumericSlot(object, "subsample", min=0, max=nrow(object))
  if(!identical(r, TRUE)) return(r)

  r = validScalarNumericSlot(object, "lts.quantile", min=0.5, max=1)
  if(!identical(r, TRUE)) return(r)

  r = validLogical(object, "verbose")
  if(!identical(r, TRUE)) return(r)

  if( (!is.character(object@calib)) ||
      (length(object@calib)!=1L) ||
      (!(object@calib)%in%c("affine","none")) ) return("Invalid slot 'calib'.")

  r = validLogical(object, "ordered")
  if(!identical(r, TRUE)) return(r)

  if(!equalOrZero(length(object@strata), nrow(object@x)))
    return("'length(strata)' must match 'nrow(x)'.")

  if(any(is.na(object@strata)))
    return("'strata' must not contain NA values.")

  if(!is.numeric(object@pstart)||length(dim(object@pstart))!=3)
    return("'pstart' must be a 3D array.")
  
  if(!all(dim(object@pstart)==c(nlevels(object@strata),
                                switch(EXPR=object@calib, affine=ncol(object@x), none=1),
                                2)))
    return("Invalid dimensions of 'pstart'.")

  if(!is.list(object@optimpar)||!identical(names(object@optimpar), optimparNames))
    return(paste("'optimpar' must be a list with elements ",
                 paste("'", compulsoryElements, "'", collapse=", ", sep=""), ".", sep=""))

  r = validScalarDoubleListElt(object@optimpar, "factr")
  if(!identical(r, TRUE)) return(r)
  
  r = validScalarDoubleListElt(object@optimpar, "pgtol")
  if(!identical(r, TRUE)) return(r)

  r = validScalarIntListElt(object@optimpar, "maxit", min=1L)
  if(!identical(r, TRUE)) return(r)

  r = validScalarIntListElt(object@optimpar, "trace", min=0L, max=6L)
  if(!identical(r, TRUE)) return(r)

  r = validScalarIntListElt(object@optimpar, "cvg.niter", min=1L)
  if(!identical(r, TRUE)) return(r)

  r = validScalarDoubleListElt(object@optimpar, "cvg.eps")
  if(!identical(r, TRUE)) return(r)

  return(TRUE)
}


##--------------------------------------------------
## validity method for 'vsn'
##--------------------------------------------------
validityVsn = function(object){
  if(any(is.na(object@coefficients))||(length(dim(object@coefficients))!=3))
    return("'coefficients' must be a 3D array and not contain NA values.")

  if(dim(object@coefficients)[3L]!=2L)
    return("'dim(coefficients)[3]' must be equal to 2.")
  
  if(length(object@sigsq)!=1L)
    return("'sigsq' must be of length 1.")

  if(length(object@hoffset)!=dim(object@coefficients)[1L])
    return("'length(hoffset)' and 'dim(coefficients)[1]' must match.")

  if(!equalOrZero(length(object@strata), length(object@mu)))
    return("'length(strata)' and 'length(mu)' must match.")

  if(length(object@strata)>0)
    if(nlevels(object@strata)!=dim(object@coefficients)[1])
      return("'nlevels(strata)' and 'dim(coefficients)[1]' must match.")

  switch(object@calib,
         affine = if(!equalOrZero(ncol(object@hx), dim(object@coefficients)[2]))
           return("'ncol(hx)' and 'dim(coefficients)[2]' must match."),
         none = if(!identical(dim(object@coefficients), c(1L, 1L, 2L)))
           return("'dim(object@coefficients)' must be 'c(1,1,2)'."),
         return("Invalid 'calib'."))

  if(!((length(object@lbfgsb)==1L)&&(is.integer(object@lbfgsb))))
    return("'lbfgsb' must be an integer of length 1.")
  
  return(TRUE)
}

##------------------------------------------------------------------------------------
## Class vsn
##------------------------------------------------------------------------------------
setClass("vsn",
  representation(
    coefficients = "array", ## 3D array: nrstrata * d * 2, with 2 parameters
                            ## for each stratum and array.
    strata = "factor", 
    mu = "numeric",
    sigsq = "numeric",
    hx = "matrix",
    lbfgsb = "integer",
    hoffset = "numeric",
    calib = "character"),
         
  prototype = list(
    coefficients = array(0, dim=c(0L, 0L, 2L)),
    strata = factor(integer(0L), levels="all"),
    mu = numeric(0L),
    sigsq = NA_real_,
    hx = matrix(0, nrow=0L, ncol=0L),
    lbfgsb = NA_integer_,
    hoffset = numeric(0L),
    datadim = c(0L, 0L),
    calib = "affine"),
         
  validity = validityVsn
)

##------------------------------------------------------------
## Class vsnInput
##------------------------------------------------------------
setClass("vsnInput",
  representation(
    x  = "matrix",     ## The n*d data matrix
                 
    reference = "vsn", ## A result from a previous fit (for reference normalization)
                 
    strata = "factor", ## Factor of length n, aligned with rows of x.
                       ## The code also recognizes a special case:
                       ## If 'strata' is of length 0, this
                       ##   is a compact way of representing the fact that there is
                       ##   only one stratum, i.e. this is equivalent to 'strata'
                       ##   of length n with all the same values.
                 
    ordered = "logical",  ## Have the rows of x already been sorted so that the
                          ## levels of 'strata' are consecutive (this is only a
                          ## non-trivial condition if there is more than one stratum.
                 
    lts.quantile = "numeric",
    subsample = "integer",
    verbose = "logical",
    calib = "character",
    pstart = "array",     ## Start parameters: see comment on slot 'coefficients'
                          ## in definition of class 'vsn'
    optimpar  = "list"),     ## See below: optimparnames
         
  prototype = list(
    x = matrix(as.numeric(NA), nrow=0L, ncol=0L),
    reference = new("vsn"),
    strata = factor(integer(0L), levels="all"),
    ordered = FALSE,
    lts.quantile = 1,
    subsample = 0L,
    verbose = TRUE,
    calib = "affine",
    pstart = array(as.numeric(NA), dim=c(1L,0L,2L)),
    optimpar = list(factr=5e7, pgtol=2e-4, 
                 maxit=60000L, trace=0L, cvg.niter=7L, cvg.eps=0)),
         
  validity = validityVsnInput)
                  

optimparNames = c("factr", "pgtol", "maxit", "trace", "cvg.niter", "cvg.eps")
