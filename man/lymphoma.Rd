\name{lymphoma}
\alias{lymphoma}
\docType{data}
\title{Intensity data for 8 cDNA chips with CLL and DLBL samples from
  the Alizadeh et al. paper in Nature 2000}
\description{8 cDNA chips from Alizadeh lymphoma paper}
\usage{data(lymphoma)}

\format{\code{lymphoma} is an exprSet containing the data from 8 chips
  out of the lymphoma data set by Alizadeh et al. (see references). Each
  chip represents two samples: on color channel 1 (CH1, Cy3, green) the
  common reference sample, and on color channel 2 (CH2, Cy5, red) the
  various disease samples. See \code{pData(lymphoma)}. The 9216x16
  matrix \code{exprs(lymphoma)} contains the background-subtracted spot
  intensities (CH1I-CH1B and CH2I-CH2B, respectively). 
}

\details{
The following chip intensity files were downloaded from the Stanford
microarray database (see link below), and following the links
\emph{Published Data}, \emph{Alizadeh AA, et al. (2000) Nature
  403(6769):503-11}
\verbatim{
-rw-r--r--    1 whuber    599246 Oct 23 13:01 lc7b019rex.DAT.gz
-rw-r--r--    1 whuber    625496 Oct 23 13:01 lc7b047rex.DAT.gz
-rw-r--r--    1 whuber    616576 Oct 23 13:01 lc7b048rex.DAT.gz
-rw-r--r--    1 whuber    611138 Oct 23 13:01 lc7b056rex.DAT.gz
-rw-r--r--    1 whuber    606123 Oct 23 13:01 lc7b057rex.DAT.gz
-rw-r--r--    1 whuber    608323 Oct 23 13:01 lc7b058rex.DAT.gz
-rw-r--r--    1 whuber    607160 Oct 23 13:01 lc7b069rex.DAT.gz
-rw-r--r--    1 whuber    598929 Oct 23 13:01 lc7b070rex.DAT.gz
}}
\references{A. Alizadeh et al., Distinct types of diffuse large B-cell lymphoma 
identified by gene expression profiling. Nature 403(6769):503-11, Feb 3, 2000.}

\source{http://genome-www5.stanford.edu/MicroArray/SMD}

\seealso{\code{\link{vsn}}}
\examples{
   data(lymphoma)
   lymphoma
   pData(lymphoma)
}

\keyword{datasets}
