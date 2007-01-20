\name{vsn2}
\alias{vsn2}
\alias{vsn2-methods}
\alias{vsn2,matrix-method}
\alias{vsn2,numeric-method}
\alias{vsn2,ExpressionSet-method}
\alias{vsnMatrix}

\title{Fit the vsn model}
\description{Fit the vsn model.
  This is the main function of this package, see also the vignette vsn.pdf.}
\usage{vsnMatrix(x,
  reference,
  strata,
  lts.quantile = 0.9,
  subsample    = as.integer(0),
  verbose      = interactive(),
  returnData   = TRUE,
  pstart,
  cvg.niter    = as.integer(7),
  cvg.eps      = 1e-3)

\S4method{vsn2}{ExpressionSet}(x, reference, strata, ...)
\S4method{vsn2}{matrix}(x, reference, strata, ...)
\S4method{vsn2}{numeric}(x, reference, strata, ...)
}
\arguments{
  \item{x}{An object contain the data to which the model is to be
    fitted. Methods exists for
    \code{\link[ExpressionSet-class]{ExpressionSet}},
    \code{matrix} and \code{numeric}.}
  \item{reference}{optional, a \code{\link[vsn-class]{vsn}} object from
    a previous fit. If this argument is specified, the data in \code{x}
    are normalized "towards" an existing set of reference arrays whose
    parameters are stored in the object \code{reference}. If this
    argument is not specified, then the data in \code{x} are normalized
    "among themselves". See details for a more formal explanation.}
  \item{strata}{optional, a factor whose length is \code{nrow(x)}. Can
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
    in the returned \code{\link[vsn-class]{vsn}} object.
    The option to set this \code{FALSE} allows saving of CPU and memory
    resources in some applications.}
  \item{pstart}{Array. optional. Can be used to specify start values
    for the iterative parameter estimation algorithm. See 
    \code{\link{vsn2trsf}} for a description of the layout of the array.}
  \item{cvg.niter}{Integer. The number of iterations to be used in the least
    trimmed sum of squares regression.}
  \item{cvg.eps}{Numeric. A convergence treshold.}
  \item{...}{Arguments that get passed on the \code{vsnMatrix}.}
}
\value{
  An object of class \code{\link[vsn-class]{vsn}}.
  The transformed data are on a glog scale to base 2. More precisely,
  the transformed data are subject to the transformation
  $\asinh(a+bx)/\log(2)-c$, where
  $\asinh(x)=\log(x+\sqrt{x^2+1})$ is also called the 'glog', and
  the constant $c$ an overall constant offset that is computed such that for
  large x the transformation approximately corresponds to the
  $\log2$. It is inconsequential for all differential expression
  calculations, but many users like to see the data on a scale that they
  are familiar with.
}

\details{
  If the \code{reference} argument is \emph{not} specified, then the model
  parameters $\mu_k$ and $\sigma$ are fit from the data in \code{x}.
  This is the mode of operation described in the 2002 Bioinformatics
  paper and that was the only option in versions 1.X of this package.
  If code{reference} is specified, the model parameters are taken from
  it. This allows for 'incremental' normalization.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")
fit = vsn2(kidney)
plot(exprs(fit), pch=".")
abline(a=0, b=1, col="red")
}

\keyword{}
