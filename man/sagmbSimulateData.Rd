\name{sagmbSimulateData}
\alias{sagmbSimulateData}
\alias{sagmbAssess}
\title{Simulate data and assess vsn's parameter estimation}
\description{Functions to validate and assess the performance of vsn
  through simulation of data.}
\usage{
sagmbSimulateData(n=8064, d=2, de=0, up=0.5, nrstrata=1)
sagmbAssess(h1, sim)
}
\arguments{
  \item{n}{Numeric. Number of probes (rows).}
  \item{d}{Numeric. Number of arrays (columns).}
  \item{de}{Numeric. Fraction of differentially expressed genes.}
  \item{up}{Numeric. Fraction of up-regulated genes among the
    differentially expressed genes.}
  \item{nrstrata}{Numeric. Number of probe strata.}
  \item{h1}{Matrix. Calibrated and transformed data, according, e.g., to
    vsn}
  \item{sim}{List. The output of a previous call to
    \code{sagmbSimulateData}, see Value}
}

\value{
  For \code{sagmbSimulateData}, a list with four components:
  \code{hy}, an \code{n x d} matrix with the true (=simulated)
  calibrated, transformed data;
  \code{y}, an \code{n x d} matrix with the simulated
  uncalibrated raw data - this is intended to be fed into
  \code{\link[vsn:vsn]{vsn}};
  \code{is.de}, a logical vector of length \code{n}, specifying
  which probes are simulated to be differentially expressed.
  \code{strata}, a factor of length \code{n}.

  For \code{sagmbSimulateData}, a number: the root mean squared
  difference between true and estimated transformed data.
}

\details{Please see the vignette.}

\references{Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann,
  Annemarie Poustka, and Martin Vingron (2003)
  "Parameter estimation for the calibration and variance stabilization
  of microarray data",
  Statistical Applications in Genetics and Molecular Biology:
  Vol. 2: No. 1, Article 3.
  http://www.bepress.com/sagmb/vol2/iss1/art3}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{ \code{\link{vsn}} }

\examples{
  sim <- sagmbSimulateData(nrstrata=4)
  ny  <- vsn(sim$y, strata=sim$strata)
  res <- sagmbAssess(exprs(ny), sim)
  res
}

\keyword{datagen}
