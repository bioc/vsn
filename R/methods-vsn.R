setMethod("predict",
  signature("vsn"),
  definition = function(object, newdata) {
    stopifnot(is.matrix(newdata))
    ## TO DO - add treatment of newdata an ExpressionSet
    ##   unfortunately the generic function does not allow dispatch on 'newdata'
    storage.mode(newdata) = "double"
    stopifnot(validObject(object))
    
    s = object@strata
    if(length(strata)==0) {
      s = rep(as.integer(1), nrow(newdata))
    } else {
      if(length(s)!=nrow(newdata))
        stop("'nrow(newdata)' must match the length of slot 'strata' in 'object'.")
    }
    
    hy = .Call("vsn2_c", newdata, as.vector(object@par), s, as.integer(2), PACKAGE="vsn")

    dimnames(hy) = dimnames(newdata)
    return(hy)
  })
       

