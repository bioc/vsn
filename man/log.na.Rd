\name{log.na}
\alias{log.na}
\title{Wrapper around the log function from base that avoids warnings
  when negative values are presented}
\description{Wrapper around the log function from base that avoids warnings
  when negative values are presented}
\usage{
log.na(x, ...)
}
\arguments{
  \item{x}{A numeric or complex vector.}
  \item{...}{Further arguments that get passed on to
    \code{\link{log}}.}
}

\value{A vector of the same length as `x' containing the transformed
  values.  For x<=0, NA is returned without warning.}

\author{Wolfgang Huber \url{http://www.dkfz.de/mga/whuber}}
\seealso{\code{\link{log}}}

\examples{
   log.na(-3:3)
}

\keyword{math}
