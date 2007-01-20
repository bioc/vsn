\name{vsn2trsf}
\alias{vsn2trsf}
\alias{predict,vsn-method}

\title{Apply the vsn transformation to data}
\description{Apply the vsn transformation to data.}
\usage{
vsn2trsf(x, p, strata)
\S4method{predict}{vsn}(object, newdata)
}

\arguments{
  \item{x, newdata}{Numeric matrix.}
  \item{object}{An object of class \code{\link[vsn-class]{vsn}} that
    contains transformation parameters and strata information.}
  \item{p}{An array with the transformation parameters. If \code{strata} is
    specified, \code{p} is expected to be a 3d array, the first
    dimension corresponds to the levels of \code{strata},
    the second dimension to the columns of \code{x} and the third dimension
    must be 2, corresponding to offsets and factors.
    If \code{strata} is missing, then the first dimension may be omitted
    and \code{p} can be a 2d array. \code{NA} values are not allowed.}
  \item{strata}{optional, a factor whose length is \code{nrow(x)}. Can
    be used for stratified normalization (i.e. separate offsets \code{a} and
    factors \code{b} for each level of \code{strata}).}
}
\value{
  A matrix of the same size as \code{x}.
}

\details{.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")
p = array(c(-0.1, -0.1, 0.003, 0.003), dim=c(1,2,2))
tx = vsn2trsf(exprs(kidney), p)
plot(tx, pch=".")
}

\keyword{}
