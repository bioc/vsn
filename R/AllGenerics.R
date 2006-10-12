## ===========================================================================
## Generic for meanSdPlot
## ---------------------------------------------------------------------------
setGeneric("meanSdPlot", function(x, ranks = TRUE, xlab = ifelse(ranks,
                                  "rank(mean)", "mean"), ylab = "sd",
                                  pch = ".", col, ...)
           standardGeneric("meanSdPlot")
)
## ===========================================================================




