
rowSds <- function(x, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<=1] = NA
  return(sqrt(rowSums(sqr(x-rowMeans(x, ...)), ...)/(n-1)))
}



## ==========================================================================
## meanSdPlot method for matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot",
            signature="matrix", definition=function(x, ranks=TRUE,
              xlab = ifelse(ranks, "rank(mean)", "mean"),
              ylab = "sd", pch  = ".", col, ...) {
                stopifnot(is.logical(ranks), length(ranks)==1, !is.na(ranks))

                ## the coloring
                pcol <- "black"
                if(!missing(col))
                  pcol <- col
                      
                if (! is(x, "matrix"))
                  stop("'x' must be a matrix or an exprSet or eSet ",
                       "(or it may inherit from these).")

                n    <- nrow(x)
                px   <- rowMeans(x, na.rm=TRUE)
                py   <- rowSds(  x, na.rm=TRUE)
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
                plot(px, py, pch=pch, xlab=xlab, ylab=ylab, col=pcol, ...)
                lines(pxl, rq.sds, col="red", type="b", pch=19)
              })
## ==========================================================================


## ==========================================================================
## meanSdPlot method for eSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot",
          signature="eSet", definition=function(x, ranks=TRUE,
            xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", col, ...) {
              if (missing(col)) {
                sel <- preproc(description(x))$vsnTrimSelection
                if (!is.null(sel)) {
                  if (!is.logical(sel) || length(sel) != nrow(exprs(x)) ||
                      any(is.na(sel)))
                    stop(paste("The element \"vsnTrimSelection\" of the preprocessing",
                               "slot of the description slot of \"x\" is not valid.",
                               "You may remove it and try again.\n"))
                  col <- ifelse(sel, "blue", "black")
                }
                meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, col=col, ...)
              }})
## ==========================================================================



## ==========================================================================
## meanSdPlot method for exprSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("meanSdPlot",
          signature="exprSet", definition=function(x, ranks=TRUE,
            xlab = ifelse(ranks, "rank(mean)", "mean"),
            ylab = "sd", pch  = ".", col, ...) {
              if (missing(col)) {
                sel <- preproc(description(x))$vsnTrimSelection
                if (!is.null(sel)) {
                  if (!is.logical(sel) || length(sel) != nrow(exprs(x)) ||
                      any(is.na(sel)))
                    stop(paste("The element \"vsnTrimSelection\" of the preprocessing",
                               "slot of the description slot of \"x\" is not valid.",
                               "You may remove it and try again.\n"))
                  col <- ifelse(sel, "blue", "black")
                }
                meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, col=col, ...)
              }})
## ==========================================================================         
