## It would be nice if we could use S4 to dispatch on the class of the second argument
##   "newdata" as well, but unfortunately the generic function "predict" from
##   package "stats" does not define a second argument. We use that function 
##   anyway to hang the following method upon, but will need to do the dispatching
##   manually.

setMethod("predict", signature("vsn"),
  function(object, newdata, strata=object@strata, log2scale=TRUE, useDataInFit=FALSE) {

    stopifnot(validObject(object), validObject(newdata),
              length(log2scale)==1L,    is.logical(log2scale),    !is.na(log2scale),
              length(useDataInFit)==1L, is.logical(useDataInFit), !is.na(useDataInFit))

    ## order matters here:
    classes = c("ExpressionSet", "NChannelSet", "RGList", "AffyBatch", "matrix", "numeric")
    for(cl in classes)
      if(is(newdata, cl))
        return(do.call(paste("predict_vsn", cl, sep="_"),
                       args = list(object=object, newdata=newdata,
                         strata=strata, log2scale=log2scale, useDataInFit)))
               
    stop(sprintf("Expecting the class of 'newdata' to be one of: %s, but found %s.",
                 paste("'", classes, "'", sep="", collapse=", "),
                 class(newdata)))
  })


##-----------------------------------------------------------------
## Here's the meat - all the other methods are just bureaucracy
##-----------------------------------------------------------------
predict_vsn_matrix = function(object, newdata, strata, log2scale, useDataInFit) {

  if(useDataInFit) {
    
    hy = exprs(object)
    
  } else {

    storage.mode(newdata) = "double"
    if(length(strata)==0L) {
      strata = rep(1L, nrow(newdata))
    } else {
      if(length(strata)!=nrow(newdata))
        stop("'nrow(newdata)' must match 'strata'.")
      strata = as.integer(int2factor(strata))
    }

    hy = vsn2trsf(x = newdata,
                  p = object@coefficients,
                  strata = strata,
                  hoffset = if(log2scale) object@hoffset else NULL,
                  calib = object@calib)
    attributes(hy) = attributes(newdata)
      
  }
  
  hy
    
}

##-----------------------------------------------------------------
predict_vsn_numeric = function(object, newdata, ...) {
  dim(newdata) = c(length(newdata), 1)
  predict_vsn_matrix(object, newdata=newdata, ...)
}

##-------------------------------------------------------------------------
predict_vsn_AffyBatch = predict_vsn_ExpressionSet = function(object, newdata, ...) {
  exprs(newdata) = predict_vsn_matrix(object, newdata=exprs(newdata), ...)
  return(newdata)
} 

##-------------------------------------------------------------------------
predict_vsn_NChannelSet = function(object, newdata, strata, log2scale, useDataInFit,
  backgroundsubtract=FALSE, foreground=c("R","G"), background=c("Rb", "Gb")) {

  if(!useDataInFit) {
    x = MatrixFromNChannelSet(newdata, backgroundsubtract=backgroundsubtract,
                              foreground=foreground, background=background)
  } else {
    x = NULL
  }

  mat = predict_vsn_matrix(object, newdata=x, strata=strata,
    log2scale=log2scale, useDataInFit=useDataInFit)

  d    = ncol(newdata)
  nrch = ncol(mat)%/%d
  stopifnot(ncol(mat)%%d==0L)
  chn = attr(object, "ChannelNames")
  stopifnot(is.character(chn), length(chn)==nrch)

  args = vector(mode="list", length=nrch+1L)
  args[[1L]] = "lockedEnvironment"
  for(i in seq_len(nrch))
      args[[i+1L]] = mat[, seq_len(d)+(i-1L)*d, drop=FALSE]
  names(args) = c("storage.mode", chn)
    
  assayData(newdata) = do.call(assayDataNew, args)
  if(!(is(newdata, "NChannelSet") && validObject(newdata)))
    stop("Failed to create a valid 'NChannelSet'.")
  return(newdata)
}

  
##-------------------------------------------------------------------------
predict_vsn_RGList = function(object, newdata, ...)
  predict_vsn_NChannelSet(object, as(newdata, "NChannelSet"), ...)
  
