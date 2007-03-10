\name{justvsn}
\alias{justvsn}
\alias{justvsn-methods}
\alias{justvsn,ExpressionSet-method}
\alias{justvsn,AffyBatch-method}
\alias{justvsn,RGList-method}

\title{Normalization with vsn}
\description{\code{justvsn} fits the vsn model to (raw microarray) data
  in an \code{\linkS4class{ExpressionSet}},
  \code{\linkS4class{AffyBatch}} or
  \code{\link{RGList}}.
  
  \code{\link{justvsn}} is a simple wrapper for
  \code{\link{vsn2}} and \code{\link{predict,vsn-method}}.
}

\usage{
justvsn(x, reference, strata, ...)
}

\arguments{
  \item{x}{An object containing the data to which the model is to be
    fitted. FIXME.}
  \item{reference, strata, ...}{As in \code{\link{vsn2}}.}
}

\value{
  FIXME.
}

\seealso{
  \code{\link{vsn2}}
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")

## use "vsn2" to produce a "vsn" object
fit = vsn2(kidney)
nkid = predict(fit, newdata=kidney)

## justvsn
nkid2 = justvsn(kidney)

stopifnot(identical(exprs(nkid), exprs(nkid2)))

}

\keyword{}
