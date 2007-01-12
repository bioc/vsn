##--------------------------------------------------
## Test whether the likelihood estimator recovers the
## correct parameters
##--------------------------------------------------

library("vsn")
options(error=recover)
dat = sagmbSimulateData(n=10000, d=1, de=0, nrstrata=1, miss=0, log2scale=TRUE)


par(mfrow=c(2,2))

## calculate the neg log likelihood:
nll = function(y, p) {
  hy = asinh(p[1]+p[2]*y)
  refh = dat$mu 
  refs = dat$sigma 
  plot(refh, hy); abline(a=0, b=1, col="blue") 
  residu = ((refh - hy)/refs)^2
  jacobi = log(p[2]/sqrt(1+(p[1]+p[2]*y)^2))
  plot(residu, jacobi)
  sum(-residu+jacobi)
}

myFun = function(p){
  cat("nll (ref):   ", nll(y=dat$y, p=p), "\n")
  cat("logLik (ref):", logLik(new("vsnInput", x=dat$y, pstart=dat$par),
                              p=p, refh=dat$mu, refsigma=dat$sigma), "\n")
}

## at the true parameters:
cat("TRUE parameters (dat$par): ", as.vector(dat$par), "\n")
myFun(as.vector(dat$par))

## fit
vr  = vsn2(dat$y, reference=new("vsn", n=length(dat$mu), refh=dat$mu, refsigma=dat$sigma),
           lts.quantile=1)

cat("\nFITTED parameters (vr@par): ", as.vector(vr@par), "\n")
myFun(as.vector(vr@par))



