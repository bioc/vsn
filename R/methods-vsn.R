#------------------------------------------------------------
# methods related to the class 'vsn'
#------------------------------------------------------------
setMethod("nrow", signature("vsn"), function(x) length(x@mu))
setMethod("ncol", signature("vsn"), function(x) switch(x@calib, affine = dim(x@coefficients)[2], none = nrow(x@hx), stop("Invalid 'calib' slot")))
setMethod("dim",  signature("vsn"), function(x) c(nrow(x), ncol(x)))

setMethod("show", signature("vsn"),
  function(object) {
    cat(class(object), sprintf("object for %d features and %d samples.\n",
      nrow(object), ncol(object)))
    if(length(object@strata)>0)
      cat(sprintf("strata: %d level%s.\n", nlevels(object@strata), c("", "s")[1+(nlevels(object@strata)>1)]))
    cat(sprintf("sigsq=%g\n", round(object@sigsq, 3)))
    if(nrow(object@hx)>0)
      cat(sprintf("hx: %d x %d matrix.\n", nrow(object@hx), ncol(object@hx)))
  })

setMethod("[", "vsn",
  function(x, i, j, ..., drop=FALSE) {
    stopifnot(missing(j), length(list(...))==0, !drop)

    x@mu = x@mu[i,drop=FALSE]
    
    if(length(x@strata)>0)
      x@strata = x@strata[i,drop=FALSE]
    if(nrow(x@hx)>0) 
      x@hx = x@hx[i,,drop=FALSE]
    
    return(x)
  })
       
setMethod("coef", signature(object="vsn"),
          function(object) object@coefficients)
setMethod("coefficients", signature(object="vsn"),
          function(object) object@coefficients)

setMethod("exprs", signature(object="vsn"),
          function(object) object@hx)

