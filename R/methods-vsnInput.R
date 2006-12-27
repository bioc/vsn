setMethod("[", "vsnInput",
  function(x, i, j, ..., drop=FALSE) {
    stopifnot(missing(j), !drop)
    x@x = x@x[i,]
    if(length(x@strata)>0)
      x@strata = x@strata[i]
    return(x)
  })
       

