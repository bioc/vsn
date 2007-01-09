setGeneric("vsn2",
  function(x, reference, strata, ...)
   standardGeneric("vsn2")) 

setGeneric("predict") ## we have imported this from package stats

setGeneric("nrow")    ## from base
setGeneric("ncol")   
setGeneric("dim")    

setGeneric("meanSdPlot",
  function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
           ylab = "sd", pch = ".", col, ...)
   standardGeneric("meanSdPlot"))




