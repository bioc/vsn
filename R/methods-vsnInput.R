setMethod("[", "vsnInput",
  function(x, i, j, ..., drop=FALSE) {
    stopifnot(!drop, length(list(...))==0)

    ## rows
    if(!missing(i)) {
      x@x = x@x[i,,drop=FALSE]
 
      if(nrow(x@reference)>0)
        x@reference = x@reference[i,,drop=FALSE]
    
      if(length(x@strata)>0)
        x@strata = x@strata[i,drop=FALSE]
    }

    ## columns
    if(!missing(j)) {
      x@x = x@x[,j,drop=FALSE]
      x@pstart = x@pstart[,j,,drop=FALSE]
    }
    
    return(x)
  })
       
setMethod("nrow", signature("vsnInput"), function(x) nrow(x@x))
setMethod("ncol", signature("vsnInput"), function(x) ncol(x@x))
setMethod("dim",  signature("vsnInput"), function(x) dim(x@x))

setMethod("show", signature("vsnInput"),
  function(object) {
    cat(class(object), sprintf("object for n=%d features and d=%d samples.\n",
      nrow(object), ncol(object)))
    if(length(object@strata)>0)
      cat(sprintf("strata: %d levels\n", nlevels(object@strata)))
    if(nrow(object@reference)==0)
      cat("No prior reference fit available.\n")
    else
      cat(sprintf("Has reference fit parameters; sigsq=%g\n", round(object@reference@sigsq, 3)))
  })


setMethod("logLik", signature(object="vsnInput"), vsnLogLik)
