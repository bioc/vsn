setGeneric("vsn2",
  def = function(x, reference, strata, ...) standardGeneric("vsn2"),
  valueClass = "vsn")

setGeneric("justvsn",
  def = function(x, reference, strata, ...) standardGeneric("justvsn"),
  valueClass = "ExpressionSet")

setGeneric("predict") ## S3 function imported from package stats
setGeneric("nrow")    ## S3 from base
setGeneric("ncol")   
setGeneric("dim")    


setGeneric("meanSdPlot",
  function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
           ylab = "sd", pch = ".", ...)
   standardGeneric("meanSdPlot"))




