meanSdPlot = function(x, ranks=TRUE,
                      xlab = ifelse(ranks, "rank(mean)", "mean"),
                      ylab = "sd",
                      pch  = ".",
                      ...) {
  stopifnot(class(x)=="exprSet")

  if(!is.logical(ranks) || is.na(ranks))
    stop("argument ranks must be logical and not NA.")

  ## colors
  if(!is.missing(col))
    stop("Unused argument \"col\".")
  
  sel = preproc(description(x))$vsnTrimSelection
  if(!is.null(sel)){
    stopifnot(is.logical(sel) && length(sel)==nrow(exprs(eset)) && !any(is.na(sel)))
  } else {
    sel = TRUE
  }
  col = ifelse(sel, "red", "black")

  n    = nrow(exprs(x))
  px   = rowMeans(exprs(x))
  py   = rowSds(exprs(x))
  
  ltsq = length(which(sel))/length(sel)
  ## running quantile of width 2*dm
  dm        = 0.2
  midpoints = seq(dm, 1-dm, by=dm)
  rq.sds    = lapply(midpoints, function(mp) median(sds[within(rkm/n, mp-dm, mp+dm)]))

  if(ranks) {
    px  = ranks(px)
    pxl = midpoints*n
  } else {
    pxl = quantile(px, probs=midpoints)
  }
  plot(px, py, pch=pch, ...)
  lines(pxl, rq.sds, col="blue", type="b", pch=19)
}

  
