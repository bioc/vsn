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
   function(x, reference, strata, ...) {
     m = exprs(x)
     fit = vsnMatrix(m,  reference, strata, optimpar=list(cvg.niter=4L),
       subsample = if(nrow(m)>30000L) 30000L else 0L, ...)
     exprs(x) = 2^fit@hx
     
      ## call RMA to do the probeset summarization,
      ## (with no background correction / normalization, that is already done)
      return(rma(x, normalize=FALSE, background=FALSE))
    })

setMethod("justvsn", "RGList",
   function(x, reference, strata, backgroundsubtract=FALSE, targets, ...) {
     if(!(is.logical(backgroundsubtract)&&(length(backgroundsubtract)==1)&&(!is.na(backgroundsubtract))))
       stop("'backgroundsubtract' must be a logical of length 1 and not NA.")
     if(!(is.matrix(x$R)&&is.matrix(x$G)&&(all(dim(x$R)==dim(x$G)))))
       stop("Invalid 'RGList' object, slots 'R' and 'G' must be matrices of same size.")
     y = cbind(x$R, x$G)
     
     if(backgroundsubtract) {
       if(!(is.matrix(x$Rb)&&is.matrix(x$Gb)&&(all(dim(x$Rb)==dim(x$Gb)))&&(all(dim(x$R)==dim(x$Rb)))))
         stop("Invalid 'RGList' object, slots 'Rb' and 'Gb' must be matrices of same size as 'R' and 'G'.")
       y = y - cbind(x$Rb, x$Gb)
     }
     
     dfc = data.frame(slide  = I(c(colnames(x$R),colnames(x$G))),
                      colour = rep(c("R", "G"), each = ncol(x)))
     row.names(dfc) = colnames(y) = paste(dfc$slide, dfc$colour, sep="-")
     dfcmeta = data.frame(labelDescription = I(c("Slide name obtained from colnames of RGList$R and RGList$G", "Label colour")),
                          row.names=c("slide", "colour"))

     fit = vsnMatrix(y, reference, strata, ...)
          
     res = new("ExpressionSet", exprs=fit@hx, phenoData = new("AnnotatedDataFrame", data=dfc, varMetadata=dfcmeta))
     return(res)
   })
