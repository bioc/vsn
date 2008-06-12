\name{vsnh}
\alias{vsnh}
\title{A function that transforms a matrix of microarray intensities.}
\description{A function that transforms a matrix of microarray
  intensities.
  This function works in conjunction with \code{\link{vsn}}.
  \code{vsn} and \code{vsnh} have been superseded by \code{\link{vsn2}} and
  the \code{\link[vsn:vsn2trsf]{predict}} method for \code{vsn} objects.
  The functions \code{vsn} and \code{vsnh} remain in the package for backward compatibility,
  but for new projects, please use \code{\link{vsn2}} and  \code{\link[vsn:vsn2trsf]{predict}}.}
}
\usage{vsnh(y, p, strata)}
\arguments{
  \item{y}{A numeric matrix containing intensity values from an array experiment. 
    It may contain NA values.}
  \item{p}{An array with the transformation parameters. If \code{strata} is
    specified, it must be a 3d array, dim(p)[1] must be greater than or
    equal to the maximum of \code{strata}, dim(p)[2] must be ncol(y),
    and dim(p)[3] must be 2. If \code{strata} is missing, then the first
    dimension may be omitted. NA values are not allowed. See Details.}
  \item{strata}{Integer vector of length nrow(y). See 
    \code{\link{vsn}} for details.}
}

\value{
  A numeric matrix of the same size as y, with the transformed data.
}
\details{
  The transformation is:
  
  vsnh(y, p, s)[k, i] = asinh( p[s[k], i, 1] + p[s[k], i, 2] * y[k, i] )
                        - log(2*p[s[1], 1, 2])
  
  where k=1:nrow(y) counts over the probes, 
  i=1:ncol(y) counts over the samples,
  p[s[k], i, 1] is the calibration offset for stratum s[k] in sample i,
  p[s[k], i, 2] is the calibration factor for stratum s[k] in sample i,
  and s[k] is the stratum of the the k-th probe.

  The constant offset - log(2*p[s[1], 1, 2]) is there to make sure that
  for large y, vsnh(y) for the first stratum on the first chip is
  approximately the same as log(y). This has no effect on
  the generalized log-ratios (glog-ratios), which are differences between
  transformed intensities, but some users are more comfortable with the
  absolute values that are obtained this way, since they are more
  comparable to the log scale.
}

\references{Variance stabilization applied to microarray data
calibration and to the quantification of differential expression,
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104.

Parameter estimation for the calibration and variance stabilization 
of microarray data, 
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, 
Annemarie Poustka, and Martin Vingron;  
Statistical Applications in Genetics and Molecular Biology (2003)
Vol. 2 No. 1, Article 3.
http://www.bepress.com/sagmb/vol2/iss1/art3.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{ \code{\link{vsn}} }

\examples{
data(kidney)
y      = exprs(kidney)
p      = array(c(-0.2, -0.1, 0.1, 0.2, 0.0026, 0.0028, 0.0030, 0.0032), dim=c(2,2,2))
strata = sample(1:2, nrow(y), replace=TRUE)
res1   = vsnh(exprs(kidney), p, strata)

res2   = asinh(p[strata,,1] + p[strata,,2] * y) - log(2*p[strata,1,2])

stopifnot(max(abs(res1 - res2)) < 1e-10)
}

