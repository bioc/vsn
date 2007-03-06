\name{vsn2}
\alias{vsn2}
\alias{vsn2-methods}
\alias{vsn2,matrix-method}
\alias{vsn2,numeric-method}
\alias{vsn2,ExpressionSet-method}
\alias{vsnMatrix}
\alias{justvsn}

\title{Fit the vsn model}
\description{\code{vsn2} fits the vsn model to (raw microarray) data
  in an \code{\linkS4class{ExpressionSet}} and
  returns a \code{\linkS4class{vsn}} object with the fit parameters and
  the transformed data matrix. \code{predict} applies the model
  and returns an \code{\linkS4class{ExpressionSet}} object.
  \code{\link{justvsn}} is a simple wrapper that 
  takes and returns either an \code{\linkS4class{ExpressionSet}} or a \code{matrix}. 
  These are the main functions of this package, see also the vignette
  \emph{Introduction to vsn}.
}
\usage{
\S4method{vsn2}{ExpressionSet}(x, reference, strata, ...)
justvsn(x, ...)

vsnMatrix(x,
  reference,
  strata,
  lts.quantile = 0.9,
  subsample    = 0L,
  verbose      = interactive(),
  returnData   = TRUE,
  pstart,
  cvg.niter    = 7L,
  cvg.eps      = 1e-3)

\S4method{vsn2}{matrix}(x, reference, strata, ...)
\S4method{vsn2}{numeric}(x, reference, strata, ...)
}
\arguments{
  \item{x}{An object contain the data to which the model is to be
    fitted. Methods exists for \code{\linkS4class{ExpressionSet}},
    \code{matrix} and \code{numeric}.}
  \item{reference}{Optional, a \code{\linkS4class{vsn}} object from
    a previous fit. If this argument is specified, the data in \code{x}
    are normalized "towards" an existing set of reference arrays whose
    parameters are stored in the object \code{reference}. If this
    argument is not specified, then the data in \code{x} are normalized
    "among themselves". See details for a more formal explanation.}
  \item{strata}{Optional, a factor whose length is \code{nrow(x)}. Can
    be used for stratified normalization (i.e. separate offsets \code{a} and
    factors \code{b} for each level of \code{strata}). If missing, all
    rows of \code{x} are assumed to come from one stratum.}
  \item{lts.quantile}{Numeric. The quantile that is used for the resistant
    least trimmed sum of squares regression. Allowed values are between
    0.5 and 1. A value of 1 corresponds to ordinary least sum of squares
    regression.}
  \item{subsample}{Integer. If specified, the model parameters are
    estimated from a subsample of the data only, the transformation is
    then applied to all data. For large datasets, this can greatly
    increase the performance at a negligible loss of precision.}
  \item{verbose}{Logical. If TRUE, some messages are printed.}
  \item{returnData}{Logical. If TRUE, the transformed data are returned
    in the returned \code{\linkS4class{vsn}} object.
    The option to set this \code{FALSE} allows saving of CPU and memory
    resources in some applications.}
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
  $\mbox{asinh}(x)=log(x+sqrt{x^2+1})$ is also called the 'glog', and
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
  If \code{reference} is specified, the model parameters are taken from
  it. This allows for 'incremental' normalization.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")

## use "vsn2" to produce a "vsn" object
fit = vsn2(kidney)
nkid = predict(fit, newdata=kidney)
plot(exprs(nkid), pch=".")
abline(a=0, b=1, col="red")

## justvsn
nkid2 = justvsn(kidney)

stopifnot(identical(exprs(nkid), exprs(nkid2)))

}

\keyword{}
