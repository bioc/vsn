\name{vsn2}
\alias{vsn2}
\alias{vsn2-methods}
\alias{vsn2,matrix-method}
\alias{vsn2,numeric-method}
\alias{vsn2,ExpressionSet-method}
\alias{vsnMatrix}

\title{Fit the vsn model}
\description{\code{vsn2} fits the vsn model to the data matrix
  in an \code{\linkS4class{ExpressionSet}} and
  returns a \code{\linkS4class{vsn}} object with the fit parameters and
  the transformed data matrix.
  The data matrix contains, typically, feature intensity readings from a
  microarray. There are also \code{vsn2} methods for numeric matrices and vectors.
  \code{\link[vsn:vsn2trsf]{predict}} applies a fitted model
  to data and returns an \code{\linkS4class{ExpressionSet}} object.
  \code{\link{justvsn}} is a simple wrapper that 
  takes and returns an \code{\linkS4class{ExpressionSet}}.
  These are the main functions of this package. An overview is given in the vignette
  \emph{Introduction to vsn}.
}
\usage{
\S4method{vsn2}{ExpressionSet}(x, reference, strata, ...)

vsnMatrix(x,
  reference,
  strata,
  lts.quantile = 0.9,
  subsample    = 0L,
  verbose      = interactive(),
  returnData   = TRUE,
  pstart,
  cvg.niter    = 5L,
  cvg.eps      = 5e-3)

\S4method{vsn2}{matrix}(x, reference, strata, ...)
\S4method{vsn2}{numeric}(x, reference, strata, ...)
}
\arguments{
  \item{x}{An object containing the data to which the model is to be
    fitted. Methods exists for \code{\linkS4class{ExpressionSet}},
    \code{matrix} and \code{numeric}.}
  \item{reference}{Optional, a \code{\linkS4class{vsn}} object from
    a previous fit. If this argument is specified, the data in \code{x}
    are normalized "towards" an existing set of reference arrays whose
    parameters are stored in the object \code{reference}. If this
    argument is not specified, then the data in \code{x} are normalized
    "among themselves". See Details for a more precise explanation.}
  \item{strata}{Optional, a factor whose length is \code{nrow(x)}. Can
    be used for stratified normalization (i.e. separate offsets \code{a} and
    factors \code{b} for each level of \code{strata}). If missing, all
    rows of \code{x} are assumed to come from one stratum.}
  \item{lts.quantile}{Numeric of length 1. The quantile that is used for the resistant
    least trimmed sum of squares regression. Allowed values are between
    0.5 and 1. A value of 1 corresponds to ordinary least sum of squares
    regression.}
  \item{subsample}{Integer of length 1. If specified, the model parameters are
    estimated from a subsample of the data of size \code{subsample}
    only, yet the fitted transformation is
    then applied to all data. For large datasets, this can substantially
    reduce the CPU time and memory consumption at a negligible loss of precision.}
  \item{verbose}{Logical. If TRUE, some messages are printed.}
  \item{returnData}{Logical. If TRUE, the transformed data are returned
    in a slot of the resulting \code{\linkS4class{vsn}} object.
    The option to set this option to \code{FALSE} allows saving of memory
    if the data are not needed.}
  \item{pstart}{Optional, array. Can be used to specify start values
    for the iterative parameter estimation algorithm. See 
    \code{\link{vsn2trsf}} for a description of the layout of the array.}
  \item{cvg.niter}{Integer. The number of iterations to be used in the least
    trimmed sum of squares regression.}
  \item{cvg.eps}{Numeric. A convergence treshold.}
  \item{...}{Arguments that get passed on to \code{vsnMatrix}.}
}
\value{
  An object of class \code{\linkS4class{vsn}}.
  The transformed data are on a glog scale to base 2. More precisely,
  the transformed data are subject to the transformation
  asinh(a+bx)/log(2)+c, where
  $\mbox{asinh}(x)/log(2)=log_2(x+sqrt{x^2+1})$ is also called the 'glog', and
  the constant c is an overall constant offset that is computed such that for
  large x the transformation approximately corresponds to the
  $\log_2$ function. The offset c is inconsequential for all differential expression
  calculations, but many users like to see the data in a range that they
  are familiar with.
}

\details{
  If the \code{reference} argument is \emph{not} specified, then the model
  parameters $\mu_k$ and $\sigma$ are fit from the data in \code{x}.
  This is the mode of operation described in the 2002 Bioinformatics
  paper and that was the only option in versions 1.X of this package.
  If \code{reference} is specified, the model parameters  $\mu_k$ and $\sigma$ are taken from
  it. This allows for 'incremental' normalization. See the vignette
  \emph{Likelihood Calculations for vsn}.
}

\seealso{\code{\link{justvsn}}, \code{\link[vsn:vsn2trsf]{predict}}}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")

fit = vsn2(kidney)                   ## fit
nkid = predict(fit, newdata=kidney)  ## apply fit

plot(exprs(nkid), pch=".")
abline(a=0, b=1, col="red")
}

\keyword{}
