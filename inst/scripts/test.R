library("vsn")
options(error=recover)

n = 5
dat = sagmbSimulateData(n=1000*(2^n), d=1, de=0, nrstrata=1, miss=0, log2scale=TRUE)
ref = new("vsn", n=length(dat$mu), refh=dat$mu, refsigma=dat$sigma)

fitpar = sapply(1:n, function(i) {
  subs = 1:(500*(2^i))
  fit = vsn2(dat$y[subs,], reference=ref[subs,], lts.quantile=1, verbose=FALSE)
  c(fit@par,
    mad(abs(dat$hy[subs,]-fit@hx)))
})

par(mfrow=c(nrow(fitpar), 1))
for(i in 1:nrow(fitpar)) {
  target = if(i<=2) dat$par[1,1,i] else 0
  plot(fitpar[i,], pch=20, ylab=letters[i], ylim=range(c(fitpar[i,], target)))
  abline(h=target, col="blue")
}

vin = new("vsnInput", x=dat$y, pstart=array(as.numeric(NA), dim=c(1,1,2)), lts.quantile = 1)

source("../../R/plotLikelihood.R")

## How does the likelihood landscape look around the true parameters?
lp = plotVsnLogLik(vin, dat$par, refh=ref@refh, refsigma=ref@refsigma, expand=1)

wm = which.max(lp$logLik)
print( lp[wm, ] )
