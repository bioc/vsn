##-----------------------------------------------------------
## mu as in the statistical model - expectation of transformed intensities 
## generate sinh(mu) according to model of Newton et al.
## shape parameter a=1, scale theta=1
##------------------------------------------------------------
sagmbSimulateData <- function(n=8064, d=2, de=0, up=0.5, nrstrata=1) {
  stopifnot(is.numeric(n),  length(n)==1, n>=1)
  stopifnot(is.numeric(d),  length(d)==1, d>=2) 
  stopifnot(is.numeric(de), length(de)==1, de>=0, de<=1)
  stopifnot(is.numeric(up), length(up)==1, up>=0, up<=1)
  stopifnot(is.numeric(nrstrata), length(nrstrata)==1)

  sigmaeps <- 0.2
  mu <- asinh(1/rgamma(n, shape=1, scale=1))

  ##------------------------------------------------------------
  ## the calibration parameters: 
  ## offsets are drawn from a uniform distribution over the interval
  ##     [-.2, +.2] * 90%-quantile of mu
  ## log(factors) from normal(0,1)
  ## overall scale from exp(runif(1)*10-5)
  ## w.l.o.g. first offset is 0 and first factor is 1
  ##------------------------------------------------------------
  Delta.a   <- 0.2 * quantile(sinh(mu), 0.9)
  pars      <- array(NA, dim=c(nrstrata, d, 2))
  pars[,,1] <- c(0, runif(d*nrstrata-1, min=-Delta.a, max=Delta.a))
  pars[,,2] <- c(1, exp(rnorm(d*nrstrata-1)))
  pars      <- pars * exp(runif(1)*10-5)

  ##------------------------------------------------------------
  ## generate simulated data
  ##------------------------------------------------------------
  is.de  <- (runif(n)<de)
  hy     <- matrix(as.numeric(NA), nrow=n, ncol=d)
  hy[,1] <- mu + rnorm(n, sd=sigmaeps)    ## array 1 is a reference
  for (j in 2:d) {
    s      <- 2 * as.numeric(runif(n)<up) - 1
    hy[,j] <- mu + as.numeric(is.de)*s*runif(n, min=0, max=2) + rnorm(n, sd=sigmaeps)
  }
  strata <- as.integer(ceiling(runif(n)*nrstrata))
  offs   <- pars[strata,,1]
  facs   <- pars[strata,,2]
  stopifnot(all(dim(facs)==dim(hy)), all(dim(offs)==dim(hy)))
  y <- offs + facs * sinh(hy)
  return(list(y=y, hy=hy, is.de=is.de, strata=strata))
}

## assess
sagmbAssess <- function(h1, sim) {
  stopifnot(all(c("y", "hy", "is.de") %in% names(sim)))
  h2    <- sim$hy
  is.de <- sim$is.de
  
  stopifnot(is.matrix(h1), is.matrix(h2), is.logical(is.de))
  n <- nrow(h1)
  d <- ncol(h1)
  if(nrow(h2)!=n || ncol(h2)!=d) 
    stop(paste("'h1' and 'h2' should be matrices of the same size, but have ",
         n, "x", d, " and ", nrow(h2), "x", ncol(h2), sep=""))
  stopifnot(length(is.de)==n)

  dh1 <- dh2 <- matrix(nrow=n, ncol=d-1)
  for (j in 2:d) {
    dh1[, j-1] <- h1[, j] - h1[, 1]
    dh2[, j-1] <- h2[, j] - h2[, 1]
  }

  nsum <- (d-1) * sum(!is.de)
  res  <- sqrt(sum((dh1[!is.de,]-dh2[!is.de,])^2) / nsum)
  return(res)
}
