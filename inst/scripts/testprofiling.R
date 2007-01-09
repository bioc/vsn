## Test concordance of maxima of profile likelihood and normal likelihood

library("vsn")
data("kidney")
x  = exprs(kidney)

nrpar = 2*ncol(x)*nrstr
istrat = as.integer(seq(0, 1, length=ncol(x)*nrstr+1)*ncol(x)*nrow(x))

invlambda  = function(y) ifelse((1:nrpar)<=(ncol(x)*nrstr), y, 
norm       = function(x) sqrt(sum(x*x))
almostequal= function(x,y) stopifnot(max(abs(x-y)/abs(x+y))<1e-10)

if(exists("vp"))
  vp = vsn2(x, lts.quantile=1)




#### CONTINUE HERE ### 



  
testfun = function(what) {
  switch(what,
    ref = {
      refh = numeric(rowMeans(x))
      refsigma = mean(diff(t(x))^2)
    },
    profile = {
      refh = refsigma = numeric(0)
    },
    stop("Bummer"))
         
  cat(what, "wait for", nrpt, "points: ")
  for (ip in 1:nrpt) {
    cat(ip, "")
    p0 = runif(nrpar)+1.2
    gro = .Call("vsn_c",      x, p0, istrat,     as.integer(1),      PACKAGE="vsn")[-1]
    grn = .Call("vsn2_point", x, p0, istrat, numeric(0), numeric(0), PACKAGE="vsn")[-1]
    ## almostequal(gro*olddinvlambdady(p0), grn*newdinvlambdady(p0))

    df[,ip,1] = gro
    df[,ip,3] = grn
    for(il in 1:nrpar) {
      dp = eps*(il==1:nrpar)
      f1o = .Call("vsn_c",      x, p0-dp, istrat, as.integer(1),          PACKAGE="vsn")[1]
      f2o = .Call("vsn_c",      x, p0+dp, istrat, as.integer(1),          PACKAGE="vsn")[1]
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

  x11(width=14, height=7)
  par(mfrow=c(2, ncol(x)*nrstr+1))
  k = c(3, 4)
  lj = list(offset=1:(ncol(x)*nrstr), factor=ncol(x)*nrstr+(1:(ncol(x)*nrstr)))
  for(j in seq(along=lj)){
    for(il in lj[[j]]) {
      plot(df[il,,k[1]], df[il,,k[2]], pch=16, xlab=paste(k[1]), ylab=paste(k[2]),
           main=sprintf("%s %s %d", what, names(lj)[j], il))
      abline(a=0, b=1, col="blue")
    }
    hist(df[lj[[j]],,k[1]]-df[lj[[j]],,k[2]], col="orange", xlab="difference",
       main=paste(names(lj)[j], "s", sep=""), breaks=30)
    abline(v=0, col="blue")
  }
}

graphics.off()
testfun("profile")
testfun("ref")
