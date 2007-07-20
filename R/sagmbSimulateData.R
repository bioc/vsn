##----------------------------------------------------------------------------
## mu as in the statistical model - expectation of transformed intensities 
## generate sinh(mu) according to model of Newton et al.
## shape parameter a=1, scale theta=1
## Since vsn 2.0, results are returned on log2 scale, rather than natural log.
##----------------------------------------------------------------------------
sagmbSimulateData <- function(n=8064, d=2, de=0, up=0.5, nrstrata=1, miss=0, log2scale=FALSE) {
  stopifnot(is.numeric(n),  length(n)==1, n>=1)
  stopifnot(is.numeric(d),  length(d)==1, d>=1) 
  stopifnot(is.numeric(de), length(de)==1, de>=0, de<=1)
  stopifnot(is.numeric(up), length(up)==1, up>=0, up<=1)
  stopifnot(is.numeric(nrstrata), length(nrstrata)==1)

  sigsq = 0.04
  mu = asinh(1/rgamma(n, shape=1, scale=1))

  ##------------------------------------------------------------
  ## the calibration parameters: 
  ##------------------------------------------------------------
  coefficients = array(runif(nrstrata*d*2L, min=-2, max=+2),
                       dim=c(nrstrata, d, 2L))

  ##------------------------------------------------------------
  ## generate simulated data:
  ##  hy = asinh(f(b)*y+a)   <=>  y = (sinh(hy)-a)/f(b)
  ##------------------------------------------------------------
  is.de  <- (runif(n)<de)
  hy     <- matrix(as.numeric(NA), nrow=n, ncol=d)
  hy[,1] <- mu + rnorm(n, sd=sqrt(sigsq))    ## array 1 is a reference
  for (j in seq_len(d)[-1]) {
    s      <- 2 * as.numeric(runif(n)<up) - 1
    hy[,j] <- mu + as.numeric(is.de)*s*runif(n, min=0, max=2) + rnorm(n, sd=sqrt(sigsq))
  }
  strata <- as.integer(ceiling(runif(n)*nrstrata))
  offs   <- coefficients[strata,,1]
  facs   <- coefficients[strata,,2]
  stopifnot(all(dim(facs)==dim(hy)), all(dim(offs)==dim(hy)))

  y = (sinh(hy)-offs)/scalingFactorTransformation(facs)
  if(miss>0)
    y[sample(length(y), length(y)*miss)] = as.numeric(NA)

  if(log2scale)
    hy = hy/log(2)
  
  return(list(y=y, hy=hy, mu=mu, sigsq=sigsq, coefficients=coefficients, is.de=is.de, strata=strata))
}

## assess
sagmbAssess <- function(h1, sim) {
  stopifnot(all(c("y", "hy", "is.de") %in% names(sim)))
  h2    <- sim$hy
  is.de <- sim$is.de

  stopifnot(is.matrix(h1), is.matrix(h2), is.logical(is.de),
            identical(is.na(h1), is.na(sim$y)), !any(is.na(h2)))
  
  n <- nrow(h1)
  d <- ncol(h1)
  if(nrow(h2)!=n || ncol(h2)!=d) 
    stop(paste("'h1' and 'h2' should be matrices of the same size, but have ",
         n, "x", d, " and ", nrow(h2), "x", ncol(h2), sep=""))
  stopifnot(length(is.de)==n)

  dh1 <- h1-rowMeans(h1, na.rm=TRUE)
  dh2 <- h2-rowMeans(h2)

  dh   <- dh1[!is.de,] - dh2[!is.de,]
  sqrt(mean(dh^2, na.rm=TRUE))
}
