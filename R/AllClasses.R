validScalar = function(ob, nm, min=0, max=+Inf) {
  s = slot(ob, nm)
  if((length(s)!=1)||any(is.na(s))||(s<min)||(s>max))
    stop(sprintf("'%s' must be a scalar with values between %g and %g.", nm, min, max))
}

validityVsnInput = function(object){
     
  validScalar(object, "subsample", min=-Inf)
  validScalar(object, "cvg.niter", min=1)
  validScalar(object, "cvg.eps")
  validScalar(object, "lts.quantile", min=0.5, max=1)
  
  if(any(is.na(object@verbose)) || (length(object@verbose)!=1))
    stop("'verbose must be of length 1 and not NA.")
  
  if(any(is.na(object@ordered)) || (length(object@ordered)!=1))
    stop("'ordered must be of length 1 and not NA.")

  if(!(length(object@strata)%in%c(0, nrow(object@x)))||any(is.na(object@strata)))
    stop("'strata' must be  factor of length 0 or nrow(x) and must not contain NAs.")
  
  if(!is.numeric(object@pstart)||length(dim(object@pstart))!=3)
    stop("'pstart' must be a 3D array.")
  
  if(!all(dim(object@pstart)==c(nlevels(object@strata), ncol(object@x), 2)))
    stop("Invalid dimensions of 'pstart'.")
  
  return(TRUE)
}


setClass("vsnInput",
  representation(
    x  = "matrix",     ## The n*d data matrix
    strata = "factor", ## Factor of length n, aligned with rows of x. Special case:
                       ##  It can also be length 0, in case there is only one stratum
    ordered = "logical",  ## Are the levels consecutive in "strata"?               
    pstart = "array",     ## Start parameters (3D array: nrstrata * d * 2)
    lts.quantile = "numeric",
    cvg.niter  = "integer",
    cvg.eps   = "numeric",
    subsample = "integer",
    verbose   = "logical"),
         
   validity = validityVsnInput)
                  
