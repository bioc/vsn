\name{vsn2trsf}
\alias{predict,vsn-method}

\title{Apply the vsn transformation to data}
\description{Apply the vsn transformation to data.}
\usage{
\S4method{predict}{vsn}(object, newdata, strata=object@strata, log2scale=TRUE)
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
  \item{strata}{Optional, a \code{factor} or \code{integer}
    that aligns with the rows of \code{newdata}; see
    the \code{strata} argument of \code{\link{vsn2}}.}
  \item{log2scale}{If \code{TRUE}, the data are returned on the glog
    scale to base 2, and an overall offset c is added (see \emph{Value}
    section of the \code{\link{vsn2}} manual page). If \code{FALSE}, the
    data are returned on the glog scale to base e, and no offset is
    added.}
}
\value{
  An \code{ExpressionSet} if \code{newdata} was an \code{ExpressionSet},
  a \code{matrix} otherwise. Dimnames and metadata are preserved.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")

## nb: for random subsampling, the 'subsample' argument of vsn
##   provides an easier way to do this
fit = vsn2(kidney[sample(nrow(kidney), 500), ])
tn = predict(fit, newdata=exprs(kidney))

}

\keyword{}
