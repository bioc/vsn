## ==========================================================================
## meanSdPlot method for matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot", signature="matrix", definition =
    function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
             ylab = "sd", pch  = ".", ...) {
      
      stopifnot(is.logical(ranks), length(ranks)==1, !is.na(ranks))

      n    <- nrow(x)
      px   <- rowMeans(x, na.rm=TRUE)
      py   <- sqrt(rowVars(x, mean=px, na.rm=TRUE))
      rpx  <- rank(px, na.last=FALSE)
                
      ## run median with centers at dm,2*dm,3*dm,... and width 2*dm
      dm        <- 0.05
      midpoints <- seq(dm, 1-dm, by=dm)
      within    <- function(x, x1, x2) { x>=x1 & x<=x2 }
      mediwind  <- function(mp) median(py[within(rpx/n, mp-dm,
                                                           mp+dm)], na.rm=TRUE)
      rq.sds    <- sapply(midpoints, mediwind)
                
      if(ranks) {
        px  <- rpx
        pxl <- midpoints*n
      } else {
        pxl <- quantile(px, probs=midpoints, na.rm=TRUE)
      }
      plot(px, py, pch=pch, xlab=xlab, ylab=ylab, ...)
      lines(pxl, rq.sds, col="red", type="b", pch=19)
})

## ==========================================================================
## meanSdPlot method for ExpressionSet and exprSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot", signature="ExpressionSet", definition =
   function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", ...)
         meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, ...))

setMethod("meanSdPlot", signature="exprSet", definition =
   function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", ...) {
     .Deprecated(msg=VSN_DEPR_MSG)
     meanSdPlot(exprs(x), ...)
   })

