vsnPlotPar = function(x, what, xlab="iter", ylab=what, ...) {
  d = ncol(exprs(x))
  switch(what,
         factors = { j = d+(1:d) },
         offsets = { j = (1:d)   },
         stop(paste("Second argument \"what\" has undefined value \", what, \"\n", sep="")))
  
  pars = preproc(description(x))$vsnParamsIter
  if(is.null(pars) || !is.matrix(pars) || nrow(pars)!=2*d)
    stop("First argument \"ny\" does not seem to be resulting from vsn.")
  
  matplot(1:ncol(pars), t(pars[j,]), type="b", pch=16, xlab=xlab, ylab=ylab, ...)
}
