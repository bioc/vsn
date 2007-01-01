setMethod("[", "vsnInput",
  function(x, i, j, ..., drop=FALSE) {
    stopifnot(missing(j), !drop)
    
    x@x = x@x[i,]
    
    if(nrow(x@reference)>0)
      x@reference = x@reference[i, ]
    
    if(length(x@strata)>0)
      x@strata = x@strata[i]
    
    return(x)
  })
       

