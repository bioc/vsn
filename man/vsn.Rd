\name{vsn}
\alias{vsn}
\title{Variance stabilization and calibration for microarray data. This is 
    the main function of this package; see also the vignette vsn.pdf.}
\description{Robust estimation of variance-stabilizing and calibrating 
    transformations for microarray data}
\usage{vsn(intensities, lts.quantile=0.5, niter=10, verbose=TRUE, pstart=NULL)}
\arguments{
  \item{intensities}{Intensity values from a microarray experiment. 
    This may be a matrix, a data frame with all numeric columns, 
    an exprSet object (from package Biobase) or an marrayRaw object 
    (from package marrayClasses).
    For the two-color printed array technology,
    each row corresponds to one spot, and the columns to the different
    wave-lengths (usually red and green, but could be any number).
    Intensities may be "background" subtracted.
    For one-color arrays, each row corresponds to a probe, and each
    column to an array. The statistical model 
    implemented by this function expects the raw intensity data, 
    not logarithmically or otherwise transformed, and not thresholded
    or "floored". NAs are not accepted.}
  \item{lts.quantile}{The quantile that is used for the resistant least
    trimmed sum of squares regression. See Details.}
  \item{niter}{The number of iterations used in the least trimmed sum of
    squares regression. See Details.}
  \item{verbose}{If TRUE, information about the progress of the
    iterative algorithm is printed.}
  \item{pstart}{Initial values for the model parameters in the iterative
    parameter estimation algorithm. If NULL, the function tries to determine
    reasonable values from the distribution of \code{intensities}.}
}

\details{
The fit algorithm assumes that after a suitable affine-linear
transformation ("calibration") of each column, for the majority of rows
it holds true that the intensities within a row can be considered
replicate measurements of each other. The parameters of the calibration
as well as those of the variance stabilizing transformation are fitted
by maximum likelihood to these data.

For a minority of rows, the proportion of which is given by
lts.quantile, this need not hold true: these are considered as outliers,
and do not contribute to the fit. The rows may correspond either to
measurement outliers, or to "differentially expressed genes".

A variance-versus-mean dependence of the intensities as described in the
reference is assumed.

Allowed values for lts.quantile are between 0.5 and 1, corresponding to
least median sum of squares regression, and to ordinary least sum of
squares regression, respectively.

\bold{Performance:} This function is slow. That is due to the nested
iteration loops of numerical optimization of the likelihood function,
and the heuristic that identifies the non-outlying data points in the
least trimmed squares regression. For large arrays with many tens of
thousands of spots, you may want to consider random subsetting: that is,
only use a subset of the e.g. 10-20,000 rows of the data matrix
\code{intensities} to fit the parameters, then apply the transformation
to all the data, using \code{\link{vsnh}}. An example for this can be
seen in the function \code{\link{normalize.AffyBatch.vsn}}, whose code
you can inspect by typing \code{\link{normalize.AffyBatch.vsn}}.
}

\value{  
An object of class \link{vsn.result}. 
}

\references{Variance stabilization applied to microarray data
calibration and to the quantification of differential expression,
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}

\seealso{\code{\link{vsn.result}}}
\examples{
data(lymphoma)

## the data in the log-log-plot
if (interactive())
  x11(width=9, height=4.5)
par(mfrow=c(1,2))

y = exprs(lymphoma)[, 3:4]
colnames(y) = pData(lymphoma)$sample[3:4]
plot(y, log="xy", pch=".", main="log-log")

## use VSN to calibrate green versus red for one array
ny <- vsn(y)
show(ny)
plot(ny@h, pch=".", main="h-h")

if (interactive())
  x11()
plot(ny)
plot(ny, what="offsets")
plot(ny, what="factors")

## this should always hold true
stopifnot(all(vsnh(y, params(ny)) == ny@h)) 
}

\keyword{robust}
