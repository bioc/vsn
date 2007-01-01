
#setGeneric("vsn", function(x)) 

setGeneric("predict") ## we have imported this from package stats
setGeneric("nrow")    ## from base

setGeneric("meanSdPlot",
  function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
           ylab = "sd", pch = ".", col, ...)
           standardGeneric("meanSdPlot")
)




