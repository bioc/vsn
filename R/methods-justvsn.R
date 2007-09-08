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
            fit = vsn2(x,  reference, strata, ...)
            exprs(x) = 2^fit@hx
            
            ## call RMA to do the probeset summarization,
            ## (with no background correction / normalization, that is already done)
            rma(x, normalize=FALSE, background=FALSE)
          })


setMethod("justvsn", "NChannelSet",
          function(x,  reference, strata, ...) {
            
            fit = vsn2(x, reference, strata, ...)
            
            d    = ncol(x)
            nrch = ncol(fit)%/%d
            stopifnot(ncol(fit)%%d==0L)
            
            args = vector(mode="list", length=nrch+1L)
            args[[1L]] = "lockedEnvironment"
            for(i in seq_len(nrch))
              args[[i+1L]] = exprs(fit)[, seq_len(d)+(i-1L)*d, drop=FALSE]
            chn = attr(fit, "ChannelNames")
            stopifnot(is.character(chn), length(chn)==nrch)
            names(args) = c("storage.mode", chn)
            
            assayData(x) = do.call("assayDataNew", args)
            return(x)
          })


setMethod("justvsn", "RGList",
          function(x, reference, strata, ...) 
            justvsn(as(x, "NChannelSet"), reference, strata, ...))
