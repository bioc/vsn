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


setMethod("vsn2", "NChannelSet",
   function(x,  reference, strata, backgroundsubtract=FALSE, 
            foreground=c("R","G"), background=c("Rb", "Gb"), ...) {
     
     y = MatrixFromNChannelSet(x, backgroundsubtract=backgroundsubtract,
                               foreground=foreground, background=background)
     res = vsnMatrix(y, reference, strata, ...)
     attr(res, "ChannelNames") = foreground
     return(res)
   })


setMethod("vsn2", "AffyBatch",
   function(x, reference, strata, subsample, ...) {
     dat = exprs(x)
     if(missing(subsample))
       subsample = if((nrow(dat)>30000L) && missing(reference)) 30000L else 0L
     vsnMatrix(dat, reference, strata, optimpar=list(cvg.niter=4L),
       subsample = subsample, ...)
    })


setMethod("vsn2", "RGList",
   function(x, reference, strata, ...)
      vsn2(as(x, "NChannelSet"), reference, strata, ...))


MatrixFromNChannelSet = function(x, backgroundsubtract, foreground, background) {
  ad = assayData(x)
  if(!all(foreground%in%channelNames(x)))
    stop("One or more elements of 'foreground' are not contained in 'channelNames(x)'.")

  ## list of matrices with the foreground values
  lmat = lapply(foreground, function(k) {
    rv = ad[[k]]
    return(rv)
  })
  ## one wide matrix with all of them next to each other
  y = do.call(cbind, lmat)
  
  if(!(is.logical(backgroundsubtract)&&(length(backgroundsubtract)==1)&&(!is.na(backgroundsubtract))))
    stop("'backgroundsubtract' must be a logical of length 1 and not NA.")
  
  if(backgroundsubtract){
    if(!all(background%in%channelNames(x)))
      stop("One or more elements of 'background' are not contained in 'channelNames(x)'.")
    stopifnot(length(background)==length(foreground))
    lmat = lapply(background, function(k) ad[[k]])
    y = y - do.call(cbind, lmat)
  }

  return(y)
}
