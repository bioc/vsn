\name{vsnh}
\alias{vsnh}
\title{A function that transforms a matrix y of microarray intensities}
\description{
  Returns h(y), where h is a calibration and variance stabilization transformation
}
\usage{
vsnh(y, p)
}
\arguments{
  \item{y}{A numeric matrix containing intensity values from an array experiment. 
    It may contain NA values.}
  \item{p}{A numeric vector with the transformation parameters, of length 2*ncol(y). 
    NA values are not allowed. See Details.}
}

\value{
  A numeric matrix of the same size as y, with the transformed data.
}
\details{
  The transformation is defined as:
  
  vsnh(y,p)[k,i] = asinh( p[i] + p[i+d] * y[k,i] )
  
  where d=ncol(y) is the number of samples, 
  i=1:d counts over the samples,
  p[1:d] are the additive calibration offsets,
  p[(d+1):(2*d)] are the calibration factors,
  and k=1:nrow(y) counts over the probes.
  
}
\references{Variance stabilization applied to microarray data
calibration and to the quantification of differential expression,
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{ \code{\link{vsn}} }

\examples{
##see vsn
}

\keyword{print}
