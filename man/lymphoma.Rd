\name{lymphoma}
\alias{lymphoma}
\docType{data}
\title{Intensity data for 8 cDNA slides with CLL and DLBL samples from
  the Alizadeh et al. paper in Nature 2000}
\description{8 cDNA chips from Alizadeh lymphoma paper}
\usage{data(lymphoma)}

\format{\code{lymphoma} is an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
  containing the data from 8 chips
  from the lymphoma data set by Alizadeh et al. (see references). Each
  chip represents two samples: on color channel 1 (CH1, Cy3, green) the
  common reference sample, and on color channel 2 (CH2, Cy5, red) the
  various disease samples. See \code{pData(lymphoma)}. The 9216x16
  matrix \code{exprs(lymphoma)} contains the background-subtracted spot
  intensities (CH1I-CH1B and CH2I-CH2B, respectively). 
}

\details{
The chip intensity files were downloaded from the Stanford
microarray database. Starting from the link below, this was done by
following the links \emph{Published Data} -> 
\emph{Alizadeh AA, et al. (2000) Nature 403(6769):503-11} -> 
\emph{Data in SMD} -> \emph{Display Data}, and selecting the following 
8 slides:
\tabular{l}{
lc7b019\cr
lc7b047\cr
lc7b048\cr
lc7b056\cr
lc7b057\cr
lc7b058\cr
lc7b069\cr
lc7b070
}
Then, the script \code{makedata.R} from the \code{scripts} subdirectory
of this package was run to generate the \R data object.
}

\references{A. Alizadeh et al., Distinct types of diffuse large B-cell 
lymphoma identified by gene expression profiling. Nature 403(6769):503-11, 
Feb 3, 2000.}

\source{http://genome-www5.stanford.edu/MicroArray/SMD}

\examples{
   data("lymphoma")
   lymphoma
   pData(lymphoma)
}

\keyword{datasets}
