vsnPlotPar = function(x, what, ylab=what, ...) {
  d = ncol(exprs(x))
  switch(what,
         factors = { j = d+(1:d) },
         offsets = { j = (1:d)   },
         stop(paste("Second argument \"what\" has undefined value \", what, \"\n", sep="")))
  
  pars = preproc(description(x))
  if(is.null(pars) || nrow(pars)!=2*d)
    stop("First argument \"ny\" does not seem to be resulting from vsn.")
  
  matPlot(1:ncol(pars), t(pars), type="b", pch=16, ylab=ylab, ...)
}
