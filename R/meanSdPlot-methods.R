## ==========================================================================
## meanSdPlot method for matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot", signature="matrix", definition =
    function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
             ylab = "sd", pch  = ".", plot = TRUE, ...) {
      
      stopifnot(is.logical(ranks), length(ranks)==1, !is.na(ranks))

      n = nrow(x)
      if(n==0L) {
        warning("In 'meanSdPlot': matrix has 0 rows, there is nothing to be done.")
        return()
      }
      
      px   = rowMeans(x, na.rm=TRUE)
      py   = sqrt(rowV(x, mean=px, na.rm=TRUE))
      rpx  = rank(px, na.last=FALSE, ties.method = "random")
                
      ## run median with centers at dm,2*dm,3*dm,... and width 2*dm
      dm        = 0.05
      midpoints = seq(dm, 1-dm, by=dm)
      within    = function(x, x1, x2) { x>=x1 & x<=x2 }
      mediwind  = function(mp) median(py[within(rpx/n, mp-dm, mp+dm)], na.rm=TRUE)
      rq.sds    = sapply(midpoints, mediwind)
                
      if(ranks) {
        res = list(rank=midpoints*n, sd=rq.sds, px=rpx, py=py)
      } else {
        res = list(quantile=quantile(px, probs=midpoints, na.rm=TRUE), sd=rq.sds, px=px, py=py)
      }
      
      if(plot) {
        plot(res$px, res$py, pch=pch, xlab=xlab, ylab=ylab, ...)
        lines(res[[1L]], res$sd, col="red", type="b", pch=19)
      }
      
      return(invisible(res))
})

## ==========================================================================
## meanSdPlot method for ExpressionSet, vsn
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot", signature="ExpressionSet", definition =
   function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", plot = TRUE, ...)
         meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, plot=plot, ...))

setMethod("meanSdPlot", signature="vsn", definition =
   function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", plot = TRUE, ...)
         meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, plot=plot, ...))
