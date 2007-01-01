setGeneric("predict") ## we have imported this from package stats

setGeneric("meanSdPlot",
  function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
           ylab = "sd", pch = ".", col, ...)
           standardGeneric("meanSdPlot")
)




