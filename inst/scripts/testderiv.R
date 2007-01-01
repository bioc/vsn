library("vsn")
data("kidney")

## Test gradient calculation in vsn C code
## nrp: number of points p0 from which to consider
## lenp: length of p0 = number of parameters = 8 = number of arrays * number of strata * 2
##
## Some of the parameter values (namely, the factors) are transformed by the
## wrapper vsn_c before being given to the C functions that calculate log-likelihood
## and its gradient. 

nrp  = 20
lenp = 8

what=c("sqr","exp")[1]

oldinvlambda = function(y) ifelse((1:8)<=4, y, log(y))
newinvlambda = function(y) ifelse((1:8)<=4, y, 
  switch(what,
         stepexp =  ifelse(y<1, log(y)+1, y),
         sqr = sqrt(y),
         exp = log(y),
         stop("Zapperlot")))

eps    <- 1e-4
istrat <- as.integer((0:4)/2*nrow(exprs(kidney)))


x  = exprs(kidney)
df = array(NA, dim=c(lenp, nrp, 4))
norm <- function(x) sqrt(sum(x*x))

almostequal=function(x,y) stopifnot(max(abs(x-y)/abs(x+y))<1e-10)

cat("Wait for", nrp, "points: ")
for (ip in 1:nrp) {
  cat(ip, "")
  p0 = runif(8)+4
  gro = .Call("vsn_c",   x, p0, istrat, as.integer(1), PACKAGE="vsn")[-1]
  grn = .Call("vsn2_point", x, p0, istrat, numeric(0), numeric(0), PACKAGE="vsn")[-1]
  #almostequal(gro*olddinvlambdady(p0), grn*newdinvlambdady(p0))

  df[,ip,1] = gro
  df[,ip,3] = grn
  for(il in 1:lenp) {
    dp = eps*(il==1:lenp)
    f1o = .Call("vsn_c",  x, p0-dp, istrat, as.integer(1), PACKAGE="vsn")[1]
    f2o = .Call("vsn_c",  x, p0+dp, istrat, as.integer(1), PACKAGE="vsn")[1]
    f1n = .Call("vsn2_point", x, p0-dp, istrat, numeric(0), numeric(0), PACKAGE="vsn")[1]
    f2n = .Call("vsn2_point", x, p0+dp, istrat, numeric(0), numeric(0), PACKAGE="vsn")[1]
    almostequal(f1o, f1n)
    almostequal(f2o, f2n)
    grn = (f2n-f1n)/norm(newinvlambda(p0+dp)-newinvlambda(p0-dp))
    gro = (f2o-f1o)/norm(oldinvlambda(p0+dp)-oldinvlambda(p0-dp))
    stopifnot(is.finite(grn), is.finite(gro))
    df[il,ip,2]  = gro
    df[il,ip,4]  = grn
  }
}
cat("\n\n")

## differences = apply(df, 1:2, sd)
## hist(differences, col="lightblue")

graphics.off(); x11(width=12, height=7);par(mfrow=c(2,4))
k = c(3, 4)
for(il in 1:lenp) {
  plot(df[il,,k[1]], df[il,,k[2]], pch=16, xlab=paste(k[1]), ylab=paste(k[2]), main=paste(il))
  abline(a=0, b=1, col="blue")
}

x11(width=8, height=4); par(mfrow=c(1,2))
for(j in list(1:4, 5:8)){
  hist(df[j,,k[1]]-df[j,,k[2]], col="orange", breaks=30)
  abline(v=0, col="blue")
}
