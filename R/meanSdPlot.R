meanSdPlot = function(x, ranks=TRUE,
                      xlab = ifelse(ranks, "rank(mean)", "mean"),
                      ylab = "sd",
                      pch  = ".",
                      col  = NULL, ...) {
  stopifnot(class(x)=="exprSet")

  if(!is.logical(ranks) || is.na(ranks))
    stop("argument ranks must be logical and not NA.")

  ## the coloring
  sel = preproc(description(x))$vsnTrimSelection
  if(!is.null(sel)){
    if(!is.logical(sel) || length(sel)!=nrow(exprs(x)) || any(is.na(sel))) 
      stop("The element \"vsnTrimSelection\" of the preprocessing slot of the description slot of \"x\" is corrupted.")
    if(!is.null(col))
      stop("Parameter \"col\" must not be specified.")
    col = ifelse(sel, "red", "black")
  } else {
    sel = TRUE
    if(is.null(col))
      col = "black"
  }

  n    = nrow(exprs(x))
  px   = rowMeans(exprs(x), na.rm=TRUE)
  py   = rowSds(  exprs(x), na.rm=TRUE)
  rpx  = rank(px, na.last=FALSE)
  
  ltsq = length(which(sel))/length(sel)
  ## running median with centers at dm, 2*dm, 3*dm, ... and width 2*dm
  dm        = 0.05
  midpoints = seq(dm, 1-dm, by=dm)
  within    = function(x, x1, x2) { x>=x1 & x<=x2 }
  rq.sds    = lapply(midpoints, function(mp) median(py[within(rpx/n, mp-dm, mp+dm)]))

  if(ranks) {
    px  = rpx
    pxl = midpoints*n
  } else {
    pxl = quantile(px, probs=midpoints)
  }
  plot(px, py, pch=pch, xlab=xlab, ylab=ylab, col=col, ...)
  lines(pxl, rq.sds, col="blue", type="b", pch=19)
}

  
