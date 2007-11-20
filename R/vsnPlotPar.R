vsnPlotPar = function(x, what, xlab="iter", ylab=what, ...) {

  stopifnot(is.character(what), length(what)==1)

  whatopts <- c("offsets", "factors")
  j <- match(what, whatopts)
  if(is.na(j))
    stop(paste("Argument 'what' must be one of ", whatopts, ", is: ", what, "\n", sep=""))

  pars = preproc(description(x))$vsnParamsIter
  if(is.null(pars) || !is.array(pars) || length(dim(pars))!=4)
    stop("Argument 'x' does not contain the vsn parameters!")

  niter <- dim(pars)[4]
  matplot(1:niter, t(matrix(pars[,,j,], ncol=niter)), type="b", pch=16, xlab=xlab, ylab=ylab, ...)
}
