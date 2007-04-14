\name{vsn2trsf}
\alias{vsn2trsf}
\alias{predict,vsn-method}

\title{Apply the vsn transformation to data}
\description{Apply the vsn transformation to data.}
\usage{
\S4method{predict}{vsn}(object, newdata, log2scale=TRUE)
vsn2trsf(x, p, strata)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{vsn}} that
    contains transformation parameters and strata information, typically
    this is the result of a previous call to \code{vsn2}.}
  \item{newdata}{Numeric \code{matrix}, \code{vector} or
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} with the data to which an
    existing fit is to be applied to.
    A typical use case is when the fit is made to a set of
    control features and then the transformation is to be applied to the
    full set of data.
    If \code{newdata} is a vector, it is considered a matrix with one
    column and \code{nrow(newdata)} rows.}
  \item{log2scale}{If \code{TRUE}, the data are returned on the glog
    scale to base 2, and an overall offset c is added (see \emph{Value}
    section of the \code{\link{vsn2}} manual page). If \code{FALSE}, the
    data are returned on the glog scale to base e, and no offset is
    added.}
  \item{x}{Numeric matrix.}
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
  For \code{vsn2trsf}, a matrix of the same size as \code{x}.
  For \code{predict}, a matrix if \code{newdata} was a matrix, and an
  \code{ExpressionSet}, \code{newdata} was an
  \code{ExpressionSet}. Dimnames and metadata are preserved.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")
param = array(c(-0.1, -0.1, 0.003, 0.003), dim=c(1,2,2))
tk = vsn2trsf(exprs(kidney), param)
plot(tk, pch=".")

## see also the 'subsample' argument of vsn for an easier way to do this:
fit = vsn2(kidney[sample(nrow(kidney), 500), ])
tn = predict(fit, newdata=exprs(kidney))

}

\keyword{}
