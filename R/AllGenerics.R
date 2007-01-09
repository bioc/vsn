setGeneric("vsn2",
  function(x, reference, strata, ...)
   standardGeneric("vsn2")) 

setGeneric("predict") ## imported from package stats
setGeneric("exprs")   ## from Biobase

setGeneric("nrow")    ## from base
setGeneric("ncol")   
setGeneric("dim")    


setGeneric("meanSdPlot",
  function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
           ylab = "sd", pch = ".", col, ...)
   standardGeneric("meanSdPlot"))




