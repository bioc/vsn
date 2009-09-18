\name{vsnPlotPar}
\alias{vsnPlotPar}
\title{Plot trajectories of calibration and transformation parameters for
  a vsn fit}
\description{Plot trajectories of calibration and transformation parameters for
  a vsn fit}
\usage{vsnPlotPar(x, what, xlab="iter", ylab=what, \dots)} 
\arguments{
  \item{x}{An object of class
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet-class}} which
    has been 
    created by the function \code{\link{vsn}}.}
  \item{what}{Character, should either be "factors" or "offsets".}
  \item{xlab}{Character, label for the x-axis.}
  \item{ylab}{Character, label for the y-axis.}
  \item{\dots}{Further arguments that get passed to plot.default.}
}
\details{The plot that is created by this function may help in assessing
  whether the parameter estimation in \code{\link{vsn}} was sufficiently
  converged.}
\value{
  The function is called for its side effect, creating a plot on the
  active graphics device.
}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link{vsn}}}
\examples{
  ## see example for vsn
}
\keyword{hplot}
