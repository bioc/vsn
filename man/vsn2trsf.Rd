\name{vsn2trsf}
\alias{predict,vsn-method}

\title{Apply the vsn transformation to data}
\description{Apply the vsn transformation to data.}
\usage{
\S4method{predict}{vsn}(object, newdata, strata=object@strata, log2scale=TRUE, useDataInFit=FALSE)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{vsn}} that
    contains transformation parameters and strata information, typically
    this is the result of a previous call to \code{vsn2}.}
  \item{newdata}{Object of class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}} (from the \code{affy} package),
    \code{\link[limma:rglist]{RGList}} (from the \code{limma} package),
    \code{matrix} or \code{numeric}, with the data to which the
    fit is to be applied to.}
  \item{strata}{Optional, a \code{factor} or \code{integer}
    that aligns with the rows of \code{newdata}; see
    the \code{strata} argument of \code{\link{vsn2}}.}
  \item{log2scale}{If \code{TRUE}, the data are returned on the glog
    scale to base 2, and an overall offset c is added (see \emph{Value}
    section of the \code{\link{vsn2}} manual page). If \code{FALSE}, the
    data are returned on the glog scale to base e, and no offset is
    added.}
  \item{useDataInFit}{If \code{TRUE}, then no transformation is
    attempted and the data stored in \code{object} is transferred
    appropriately into resulting object, which otherwise preserves the
    class and metadata of \code{newdata}. This option exists to increase
    performance in constructs like
\preformatted{
       fit = vsn2(x, ...)
       nx = predict(fit, newdata=x)
  }
    and is used, for example, in the \code{\link{justvsn}} function.
  }
}

\value{
  An object typically of the same class as \code{newdata}. There are two
  exceptions: if \code{newdata} is an
  \code{\link[limma:rglist]{RGList}}, the return value is an
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, and
  if \code{newdata} is numeric, the return value is a \code{matrix} with 1
  column.
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
