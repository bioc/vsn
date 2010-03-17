\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method with expresso}
\description{Wrapper for \code{\link{vsn2}} to be used as a normalization
  method with the expresso function of the package affy. The expresso
  function is deprecated, consider using \code{\link{justvsn}}
  instead. The normalize.AffyBatch.vsn can still be useful on its own,
  as it provides some additional control of the normalization process
  (fitting on subsets, alternate transform parameters). }
\usage{
  normalize.AffyBatch.vsn(
     abatch,
     reference,
     strata = NULL,
     subsample = if (nrow(exprs(abatch))>30000L) 30000L else 0L,
     subset,
     log2scale = TRUE,
     log2asymp=FALSE,
     ...)}
\arguments{
  \item{abatch}{An object of type
    \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{reference}{Optional, a 'vsn' object from a previous fit. If this
    argument is specified, the data in 'x' are normalized "towards" an
    existing set of reference arrays whose parameters are stored in the
    object 'reference'. If this argument is not specified, then the data
    in 'x' are normalized "among themselves". See  \code{\link{vsn2}}
    for details.}
  \item{strata}{The 'strata' functionality is not supported, the
    parameter is ignored.}
  \item{subsample}{Is passed on to \code{\link{vsn2}}.}
  \item{subset}{This allows the specification of a subset of expression
    measurements to be used for the vsn fit. The transformation with the
    parameters of this fit is then, however, applied to the whole
    dataset. This is useful for excluding expression measurements that
    are known to be differentially expressed or control probes that may
    not match the vsn model, thus avoiding that they influence the
    normalization process. This operates at the level of probesets, not
    probes. Both 'subset' and 'subsample' can be used together.}
  \item{log2scale}{If TRUE, this will perform a global affine transform on
    the data to put them on a similar scale as the original
    non-transformed data. Many users prefer this.  Fold-change estimates
    are not affected by this transform. In some situations,
    however, it may be helpful to turn this off, e.g., when comparing
    independently normalized subsets of the data.}
  \item{log2asymp}{If TRUE, this will perform a global affine transform on
    the data to make the generalized log (asinh) transform be
    asymptotically identical to a log base 2 transform. Some people find
    this helpful. Only \bold{one} of 'log2scale' or 'log2asymp' can be
    set to TRUE. Fold-change estimates are not affected by this
    transform.}
  \item{...}{Further parameters for \code{\link{vsn2}}.}
}
 
\details{Please refer to the \emph{Details} and \emph{References}
  sections of the
  man page for \code{\link{vsn2}} for more details about this method.

  \bold{Important note}: after calling \code{\link{vsn2}}, the function
  \code{normalize.AffyBatch.vsn} \bold{exponentiates} the data (base 2).
  This is done in order to make the behavior of this function similar to the
  other normalization methods in affy. That packages uses the convention
  of taking the logarithm to base in
  subsequent analysis steps (e.g. in \code{\link[stats]{medpolish}}).
}

\value{An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.
  The \code{vsn} object returned, which can be used as \code{reference} for
  subsequent fits, is provided by
  \code{description(abatch)@preprocessing$vsnReference}.
}
\author{D. P. Kreil \url{http://bioinf.boku.ac.at/},
        Wolfgang Huber}

\seealso{\code{\link{vsn2}}}
\examples{
## Please see vignette.
}

