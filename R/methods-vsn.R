#------------------------------------------------------------
# methods related to the class 'vsn'
#  (see below for methods for the generic function 'vsn2')
#------------------------------------------------------------
setMethod("predict", signature("vsn"),
  function(object, newdata) {

    stopifnot(validObject(object))

    ## Treat the different types of newdata. We do this manually here,
    ## unfortunately the generic function stats::predict does not allow
    ## dispatch on 'newdata' via S4
    if(is.vector("newdata"))
      dim(newdata)=c(length(newdata), 1)

    if(is(x, "ExpressionSet")){
      stopifnot("'newdata' of class 'ExpressionSet' is not yet implemented, please supply a matrix.")
    }
    
    stopifnot(is.matrix(newdata))
    if(storage.mode(newdata) != "double")
      storage.mode(newdata) = "double"
    
    s = object@strata
    if(length(s)==0) {
      s = rep(as.integer(1), nrow(newdata))
    } else {
      if(length(s)!=nrow(newdata))
        stop("'nrow(newdata)' must match the length of slot 'strata' in 'object'.")
      s = as.integer(s)
    }
    
    hy = .Call("vsn2_trsf", newdata, as.vector(object@par), s, PACKAGE="vsn")

    dimnames(hy) = dimnames(newdata)
    return(hy)
  })
       

setMethod("nrow", signature("vsn"), function(x) x@n)
setMethod("ncol", signature("vsn"), function(x) dim(x@par)[2])
setMethod("dim",  signature("vsn"), function(x) c(x@n, dim(x@par)[2]))

setMethod("show", signature("vsn"),
  function(object) {
    cat(class(object), sprintf("object for n=%d features and d=%d samples.\n",
      nrow(object), ncol(object)))
    if(length(object@strata)>0)
      cat(sprintf("strata: %d level%s.\n", nlevels(object@strata), c("", "s",)[1+(nlevels(object@strata)>1)]))
    if(nrow(object@data)>0)
      cat(sprintf("data: %d x %d matrix.\n", nrow(object@data), ncol(object@data)))
    if(length(object@refh)>0)
      cat(sprintf("refsigma: %g.\n", round(object@refsigma, 3)))
  })

setMethod("[", "vsn",
  function(x, i, j, ..., drop=FALSE) {
    stopifnot(missing(j), length(list(...))==0, !drop)

    if(length(x@strata)>0)
      x@strata = x@strata[i,drop=FALSE]
    if(length(x@refh)>0)
      x@refh = x@refh[i,drop=FALSE]
    if(nrow(x@data)>0)
      x@data = x@data[i,,drop=FALSE]

    return(x)
  })
       
setMethod("exprs", signature(object="vsn"),
          function(object) object@data)

#------------------------------------------------------------
# methods for the generic function 'vsn2'
#------------------------------------------------------------
setMethod("vsn2", "matrix", vsnMatrix)

setMethod("vsn2", "numeric",
   function(x, reference, strata, ...)
      vsnMatrix(as.matrix(x, ncol=1), reference, strata, ...))

setMethod("vsn2", "ExpressionSet",
   function(x, reference, strata, ...)
      vsnMatrix(exprs(x), reference, strata, ...))