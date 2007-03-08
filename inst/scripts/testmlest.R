##------------------------------------------------------
## Test whether the likelihood estimator recovers the
## correct parameters, and plot the likelihood landscape
##--------------------------------------------------

library("vsn")
options(error=recover)

## Generate data
n = 500*(2^(1:7))
set.seed(569)
dat = sagmbSimulateData(n=n[length(n)], d=1, de=0, nrstrata=1, miss=0, log2scale=TRUE)
ref = new("vsn", n=length(dat$mu), refh=dat$mu, refsigma=dat$sigma, par=dat$par)
vin = new("vsnInput", x=dat$y, pstart=array(as.numeric(NA), dim=c(1,1,2)), lts.quantile = 1)

## series of fits with ascending amounts of data
fitpar = sapply(n, function(k) {
  subs = 1:k
  fit = vsn2(dat$y[subs,], reference=ref[subs,], lts.quantile=1, verbose=FALSE)
  m = mad(abs(dat$hy[subs,]-fit@hx))
##  browser()
  c(fit@par,m)
})

graphics.off()

## do the estimated parameters converge to the true ones when the data is increased?
x11(height=8, width=5)
par(mfrow=c(nrow(fitpar), 1))
for(i in 1:nrow(fitpar)) {
  target = if(i<=2) dat$par[1,1,i] else 0
  plot(sqrt(n), fitpar[i,], pch=20, xlab=expression(sqrt(n)), main="convergence", 
       ylab=letters[i], ylim=range(c(fitpar[i,], target)))
  abline(h=target, col="blue")
}

## pedestrian calculation of logLik
myLLfun = function(object, p, refh, refsigma) {
  res = matrix(as.numeric(NA), nrow=1, ncol=ncol(p))
  stopifnot(nrow(p)==2)
  for(j in 1:ncol(p))  {
    hy = asinh(p[1,j]+p[2,j]*object@x)
    residu = 0.5*((refh - hy)/refsigma)^2
    jacobi = log(p[2,j]/sqrt(1+(p[1,j]+p[2,j]*object@x)^2))
    res[, j] = sum(-residu+jacobi)
  }
  return(res)
}

## How does the likelihood landscape look around the true parameters?
ex = c(0.05, 0.05)
x11(xpos=400)
lp1 = plotVsnLogLik(vin, dat$par, refh=ref@refh, refsigma=ref@refsigma, expand=ex)

x11(xpos=800)
lp2 = plotVsnLogLik(vin, dat$par, refh=ref@refh, refsigma=ref@refsigma, expand=ex, fun=myLLfun)
if(max(abs(lp1$logLik-lp2$logLik))>1e-7) warning("lp1 and lp2 are different!")

wm = which.max(lp1$logLik)
cat(sprintf("True parameters   %11g %11g\n", dat$par[1], dat$par[2]))
cat(sprintf("Maximum of logLik %11g %11g\n\n", lp1[wm,1], lp1[wm,2]))
