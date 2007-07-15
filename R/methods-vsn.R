#------------------------------------------------------------
# methods related to the class 'vsn'
#  (see below for methods for the generic function 'vsn2')
#------------------------------------------------------------
setMethod("predict", signature("vsn"),
  function(object, newdata, strata=object@strata, log2scale=TRUE) {

    stopifnot(validObject(object))

    ## Treat the different types of newdata. We do this manually here,
    ## unfortunately the generic function stats::predict does not allow
    ## to dispatch on 'newdata' via S4
    if(is(newdata, "ExpressionSet")){
      res = newdata
      exprs(res) = predict(object, exprs(newdata), strata=strata, log2scale=log2scale)
      return(res)
    }
    
    ## If we get here, 'newdata' must be either a matrix or a vector
    if(is.vector(newdata))
      dim(newdata)=c(length(newdata), 1)

    returnClass = "matrix"
    stopifnot(is.matrix(newdata))
    if(storage.mode(newdata) != "double")
      storage.mode(newdata) = "double"
    
    if(length(strata)==0) {
      strata = rep(as.integer(1), nrow(newdata))
    } else {
      if(length(strata)!=nrow(newdata))
        stop("'nrow(newdata)' must match 'strata'.")
      strata = as.integer(int2factor(strata))
    }
    
    hy = .Call("vsn2_trsf", newdata, as.vector(object@coefficients), strata, PACKAGE="vsn")

    if(log2scale)
      hy = trsf2log2scale(hy, object@coefficients)
    
    dim(hy) = dim(newdata)
    dimnames(hy) = dimnames(newdata)
    return(hy)
  })
       

setMethod("nrow", signature("vsn"), function(x) length(x@mu))
setMethod("ncol", signature("vsn"), function(x) dim(x@coefficients)[2])
setMethod("dim",  signature("vsn"), function(x) c(nrow(x), ncol(x)))

setMethod("show", signature("vsn"),
  function(object) {
    cat(class(object), sprintf("object for n=%d features and d=%d samples.\n",
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

##------------------------------------------------------------
## methods for the generic function 'vsn2'
##------------------------------------------------------------
setMethod("vsn2", "matrix", vsnMatrix)

setMethod("vsn2", "numeric",
   function(x,  reference, strata, ...)
      vsnMatrix(as.matrix(x, ncol=1), reference, strata, ...))

setMethod("vsn2", "ExpressionSet",
   function(x,  reference, strata, ...)
      vsnMatrix(exprs(x),  reference, strata, ...))
