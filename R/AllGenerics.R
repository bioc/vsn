setGeneric("vsn2",
  function(x, reference, strata, ...)
   standardGeneric("vsn2")) 

setGeneric("justvsn",
  function(x, reference, strata, ...)
   standardGeneric("justvsn")) 

setGeneric("predict") ## S3 function imported from package stats
setGeneric("nrow")    ## S3 from base
setGeneric("ncol")   
setGeneric("dim")    


setGeneric("meanSdPlot",
  function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
           ylab = "sd", pch = ".", ...)
   standardGeneric("meanSdPlot"))




