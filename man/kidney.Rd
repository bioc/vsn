\name{kidney}
\alias{kidney}
\docType{data}
\title{Intensity data for 1 cDNA slide with two adjacent tissue samples
  from a nephrectomy (kidney)}
\description{Intensity data for 1 cDNA slide with two adjacent tissue samples
  from a nephrectomy (kidney)}
\usage{data(kidney)}

\format{\code{kidney} is an exprSet containing the data from 1 cDNA
  chip. The 8704x2 matrix \code{exprs(kidney)} contains the
  spot intensities for the red (635 nm) and green color channels
  (532 nm) respectively. For each spot, a background estimate from a
  surrounding region was subtracted.
}

\details{The chip was produced in 2001 by Holger Sueltmann at
  the Division of Molecular Genome Analysis at the German Cancer
  Research Center in Heidelberg.}

\references{Huber W, Boer JM, von Heydebreck A, Gunawan B, Vingron M,
  Fuzesi L, Poustka A, Sueltmann H. Transcription profiling of renal
  cell carcinoma. Verh Dtsch Ges Pathol. 2002;86:153-64.
  PMID: 12647365}

\seealso{\code{\link{vsn}}}
\examples{
 data(kidney)
 plot(exprs(kidney), pch=".", log="xy")
 abline(a=0,b=1,col="blue")  
}
\keyword{datasets}
