\name{getIntensityMatrix}
\alias{getIntensityMatrix}
\title{Extract a matrix of microarray intensities from various kinds of objects}
\description{Extract a matrix of microarray intensities from various
  kinds of objects. This function is called by the function
  \code{\link[vsn:vsn]{vsn}}, and is normally not called by the user.}
\usage{getIntensityMatrix (intensities, verbose)}
\arguments{
\item{intensities}{Object of class
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
  with raw intensity values from a microarray experiment. Alternatively,
  this may be a \code{\link{matrix}}, a \code{\link{data.frame}} with
  all numeric columns or an
  \code{\link[marray:marrayRaw-class]{marrayRaw}} object.}
\item{verbose}{Logical. If TRUE, some messages are printed.}
}
\value{
  A numeric matrix.
}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link[vsn:vsn]{vsn}}}
\examples{
  ## see the vignette
}
\keyword{manip}
