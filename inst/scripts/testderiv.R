library(vsn)
data(kidney)

eps    <- 1e-4
istrat <- as.integer((0:4)/2*nrow(exprs(kidney)))

norm <- function(x) sqrt(sum(x*x))

nrp  <- 30
lenp <- 8
facs <- 5:8
df   <- array(NA, dim=c(2, lenp, nrp))

for (ip in 1:nrp) {
  p0 <- c(runif(4)*20-10, runif(4)*6-3)
  for(il in 1:lenp) {
    dp     <- rep(0, lenp)
    dp[il] <- eps
    p      <- p0
    p[facs] <- exp(p[facs]) + eps
    f1   <- .Call("vsnc", exprs(kidney), p-dp, istrat, FALSE)[1]
    f2   <- .Call("vsnc", exprs(kidney), p+dp, istrat, FALSE)[1]
    gr   <- .Call("vsnc", exprs(kidney), p,    istrat, FALSE)[-1]
    stopifnot(!any(is.na(c(f1, f2, gr))))

    gr[facs] <- gr[facs]/p[facs] ## chain rule
    df[1, il, ip]  <- sum(gr*dp)/norm(dp)
    df[2, il, ip]  <- (f2[1]-f1[1])/(2*norm(dp))
    cat(df[1, il, ip], df[2, il, ip], df[1, il, ip]/df[2, il, ip], "\n", sep="\t")
  }
}
par(mfrow=c(2,4))
for(il in 1:lenp) {
  plot(df[1, il, ], df[2, il, ], pch=16, main=paste(il))
  abline(a=0, b=1, col="blue")
}