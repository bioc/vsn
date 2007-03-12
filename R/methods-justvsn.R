#------------------------------------------------------------
# methods for the generic function 'justvsn'
#------------------------------------------------------------

setMethod("justvsn", "ExpressionSet",
   function(x,  reference, strata, ...) {
      fit = vsnMatrix(exprs(x), reference, strata, ...)
      exprs(x) = fit@hx
      return(x)
    })

setMethod("justvsn", "AffyBatch",
   function(x,  reference, strata, ...) {
      fit = vsnMatrix(exprs(x),  reference, strata, ...)

      ## call RMA to do the probeset summarization,
      ## (with no background correction / normalization, that is already done)

      warning("Not Implemented")
      
      return(x)
    })

setMethod("justvsn", "RGList",
   function(x,  reference, strata, backgroundsubtract=FALSE, ...) {
     if(!(is.logical(backgroundsubtract)&&(length(backgroundsubtract)==1)&&(!is.na(backgroundsubtract))))
       stop("'backgroundsubtract' must be a logical of length 1 and not NA.")
     
      # fit = vsnMatrix(exprs(x),  reference, strata, ...)

      warning("Not Implemented")
      
      return(x)
    })
