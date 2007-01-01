setMethod("predict",
  signature("vsn"),
  definition = function(object, newdata) {

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
    
    hy = .Call("vsn2_c", newdata, as.vector(object@par), s, as.integer(2), PACKAGE="vsn")

    dimnames(hy) = dimnames(newdata)
    return(hy)
  })
       

setMethod("nrow",
  signature("vsn"),
  definition = function(x) x@n)

setMethod("show",
  signature("vsn"),
  definition = function(object) {
    cat(class(object), sprintf("object for n=%d features and d=%d samples.\n",
      object@n, dim(object@par)[2]))
    if(length(object@strata)>0)
      cat(sprintf("strata: %d levels\n", nlevels(object@strata)))
    if(nrow(object@data)>0)
      cat(sprintf("data: %d x %d matrix\n", nrow(object@data), ncol(object(data))))
    if(length(object@refh)>0)
      cat(sprintf("refsd: %g\n", object@refsd))
  })

