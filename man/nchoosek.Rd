\name{nchoosek}
\alias{nchoosek}
\title{List all subsets of size k from n objects}
\description{List all subsets of size k from n objects}
\usage{nchoosek(n, k)}
\arguments{
  \item{n}{Integer}
  \item{k}{Integer}
}

\value{  
A matrix with k rows and as many columns as there are distinct
subsets. Each column contains integers between 1 and n, and represents a
subset. }

\seealso{\code{\link{choose}}}
\examples{
nchoosek(5,3)
}

\keyword{math}
