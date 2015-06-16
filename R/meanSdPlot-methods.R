## ==========================================================================
## meanSdPlot method for matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot", signature="matrix", definition =
    function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
             ylab = "sd", pch, plot = TRUE, bins = 50, ...) {
      
      stopifnot(is.logical(ranks), length(ranks)==1, !is.na(ranks))

      n = nrow(x)
      if (n == 0L) {
        warning("In 'meanSdPlot': matrix has 0 rows, there is nothing to be done.")
        return()
      }
      if (!missing(pch)) {
        warning("In 'meanSdPlot': 'pch' is ignored.")
      }
      
      px   = rowMeans(x, na.rm=TRUE)
      py   = sqrt(rowV(x, mean=px, na.rm=TRUE))
      rpx  = rank(px, na.last=FALSE, ties.method = "random")
                
      ## run median with centers at dm, 2*dm, 3*dm,... and width 2*dm
      dm        = 0.025
      midpoints = seq(dm, 1-dm, by = dm)
      within    = function(x, x1, x2) { x>=x1 & x<=x2 }
      mediwind  = function(mp) median(py[within(rpx/n, mp-2*dm, mp+2*dm)], na.rm=TRUE)
      rq.sds    = sapply(midpoints, mediwind)
                
      res = if(ranks) {
        list(rank=midpoints*n, sd=rq.sds, px=rpx, py=py)
      } else {
        list(quantile=quantile(px, probs=midpoints, na.rm=TRUE), sd=rq.sds, px=px, py=py)
      }
      
      ## plot(res$px, res$py, pch=pch, xlab=xlab, ylab=ylab, ...)
      ## smoothScatter(res$px, res$py, xlab=xlab, ylab=ylab, ...)
      ## lines(res[[1L]], res$sd, col="red", type="b", pch=19)
      fmt = function() {
      function(x) format(round(x, 0), nsmall = 0L, scientific = FALSE)
      }
      res$gg = ggplot(data.frame(px = res$px, py = res$py),
            aes(x = px, y = py)) + xlab(xlab) + ylab(ylab) +
            stat_binhex(bins = bins, ...) +
            scale_fill_gradient(name = "count", trans = "log", labels = fmt()) + 
            geom_line(aes(x = x, y = y), data = data.frame(x = res[[1]], y = res$sd), color = "red")
              
      if (plot) print(res$gg)
      
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

setMethod("meanSdPlot", signature="MAList", definition =
   function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", plot = TRUE, ...)
         meanSdPlot(with(x, cbind(A+M/2, A-M/2)),
                              ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, plot=plot, ...))
