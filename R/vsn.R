##-----------------------------------------------------------------
## Robust calibration and variance stabilization
## (C) Wolfgang Huber 2002
## w.huber@dkfz.de
##-----------------------------------------------------------------
require(Biobase) || stop("can't load without package \"Biobase\"")

##------------------------------------------------------------
## vsn: the main function of this library
##------------------------------------------------------------
vsn <-  function(intensities, lts.quantile=0.75, niter=10, verbose=TRUE, pstart=NULL) {
  ## check lts.quantile for plausibility
  if (!is.numeric(lts.quantile) || (length(lts.quantile)!=1) || (lts.quantile<=0.5) || lts.quantile>1)
    stop(paste("invalid argument lts.quantile, expecting scalar between 0.5 and 1, but found", lts.quantile))
  ## check niter for plausibility
  if (!is.numeric(niter) || (length(niter)!=1) || (niter<4))
    stop(paste("invalid argument niter, expecting a single positive integer >=4 but found", niter))
  ## check pstart for plausibility
  if (!is.null(pstart))
    if (!is.numeric(pstart) || length(pstart) != 2*d || any(is.na(pstart)) ||
        any(pstart[(ncol(intensities)+1):(2*ncol(intensities))]<=0))
      stop(paste("invalid argument pstart, expecting numeric vector of length", 2*ncol(intensities), "but found", pstart))
  if (!is.logical(verbose))
    stop("argument verbose must be of type logical.")
  
  ## coerce different possible objects that carry the intensity data into a matrix
  y = switch(class(intensities),
     matrix     = {  if (!is.numeric(intensities))
                        stop("argument intensities must be numeric")
                     intensities
                   },
    data.frame  = {  if (!all(sapply(intensities, is.numeric)))
                       stop("argument intensities has non-numeric columns")
                     as.matrix(intensities)
                  },
    exprSet     = { tmp <- exprs(intensities)
                    colnames(tmp) <- rownames(pData(intensities))
                    tmp
                  },
    marrayRaw   = { nrslides <- ncol(intensities@maRf)
                    nrspots  <- nrow(intensities@maRf)
                    if (verbose)
                      cat(sprintf("Converting marrayRaw (%d spots, %d slides) to %dx%d matrix.\n",
                                  as.integer(nrspots), as.integer(nrslides),
                                  as.integer(nrspots), as.integer(2*nrslides)),
                          "Gf-Gb in odd columns, Rf-Rb in even columns.\n")
                    tmp = matrix(NA, nrow=nrspots, ncol=2*nrslides)
                    tmp[, (1:nrslides)*2-1 ] <- intensities@maGf - intensities@maGb
                    tmp[, (1:nrslides)*2   ] <- intensities@maRf - intensities@maRb
                    tmp
                   },
    stop(paste("argument intensities has class ", class(intensities),
               ". Permitted are: matrix, data.frame, exprSet, marrayRaw", sep=""))
    )
  
  if (any(is.na(y)) && !na.rm)
    stop("argument intensities must not contain NAs.\n")

  if (verbose)
    cat(sprintf("Called vsn on %dx%d matrix, with lts.quantile=%4.2f and niter=%d. This may take a while.\n",
                 nrow(y), ncol(y), lts.quantile, as.integer(niter)));
  
  ly <- asly <- res <- array(NA, dim=dim(y))
  nrs <- ncs <- as.integer(0)   
  ssq    <- 0.0
  last.p <- numeric(2*ncol(y))
  
  #----------------------------------------------------------
  # affine linear transformation according to p
  # output: transformed data, same size as y
  #----------------------------------------------------------
  at <- function(p) {
    stopifnot(all(!is.na(p)) & 2*ncol(y)==length(p))   # be paranoid
    offs <- matrix(p[         1  :  ncol(y) ], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
    facs <- matrix(p[(ncol(y)+1):(2*ncol(y))], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
    return(asinh(offs + facs * y))
  }
  
  #----------------------------------------------------------
  # Profile log likelihood of the model
  #
  #    asinh(a_i + b_i * y_ki) = m_k + epsilon_ki
  #
  # where k=1..n, i=1..d, y is a n-by-d matrix,
  # (a_i, b_i) are real parameters to be estimated
  # epsilon_ki are i.i.d N(0, sigma^2)
  # and sigma, and all m_k are estimated from the data
  # ("profiling")
  #
  # Argments:
  # p numeric vector of length 2*d, with elements a_1,...,a_d, b_1,...,b_d
  #
  # Variables used from the enclosing environment:
  # sy     the matrix of the y_ki
  # ly     the matrix of the a_i + b_i * y_ki
  # asly   the matrix of the asinh(a_i + b_i * y_ki) 
  # res    the matrix of the residuals epsilon_ki
  # nrs, ncs  number of rows and columns of the above matrices
  # ssq    sum of squared residuals sum(res^2)
  # last.p a place where llat() stores its call arguments. From this, when grllat() 
  #        is called it can double-check whether it is indeed also called 
  #        with the same arguments (see below for details)
  # 
  # Return value:
  # Profile Log-Likelihood
  #----------------------------------------------------------
  llat <- function(p) {
    ##if(any(is.na(p)) || any(is.na(sy)) || 2*ncol(sy)!=length(p)) {
    ##	cat("Sanity checks for parameters p failed at entry of function llat()!\n")
    ##	cat("Starting browser() to enable debugging!\n")
    ##    browser()
    ## }
    
    ## make sure the variables that are going to be assigned by "<<-" actually exist.
    stopifnot(all(unlist(lapply(c("ly","asly","res","nrs", "ncs", "ssq","last.p"),
                                exists, envir=parent.env(environment())))))
    
    nrs  <<- nrow(sy)    ## these asignments will be made in the enclosing
    ncs  <<- ncol(sy)    ## environment, i.e. the one of vsn()
    offs <- p[    1  :  ncs  ]
    facs <- p[(ncs+1):(2*ncs)]
    ly   <<- offs + facs * t(sy)    ## sy is an nxd matrix, offs and facs are d-vectors
                                    ## the cycling goes column-wise, hence transpose!
    
    asly <<- t(asinh(ly))
    last.p <<- p
    
    ## residuals
    res <<- asly - rowMeans(asly)
    ssq <<- sum(res*res)
    return( nrs*ncs/2*log(ssq) - sum(log(facs/sqrt(1+ly*ly))) )
  }
  #--------------------------------------------------------------
  # Gradient of the profile log likelihood
  #--------------------------------------------------------------
  grllat <- function(p) {
    ## Generally, optim() will call the gradient of the objective function (gr)
    ## immediately after a call to the objective (fn) at the same parameter values.
    ## Anyway, we like to doublecheck
    if(any(p!=last.p)) {
	cat("Sanity checks for parameters p failed at entry of function grllat()!\n")
	cat("Starting browser() to enable debugging!\n")
        browser()
    }
    dhda      <-   1/sqrt(1+ly*ly)
    dlndhdyda <-   -ly/(1+ly*ly)
    gra       <-   nrs*ncs/ssq*res*t(dhda) - t(dlndhdyda)
    return(c(
      colSums(gra), 
      colSums(gra*sy) - nrs/p[(ncs+1):(ncs*2)]))
  }

  #---------------------------------------------------------------------
  # guess a parameter scale, and set boundaries for optimization,
  # and if they are not supplied as arguments, set the start parameters
  #---------------------------------------------------------------------
  pscale  <- plower <- numeric(2*ncol(y))
  pscale[1:ncol(y)] <- 1
  for (j in 1:ncol(y)) {
    pscale[ncol(y)+j] <- 1/diff(quantile(y[,j], probs=c(0.25, 0.75)))
  }
  # lower boundary for 'factors': a small positive value
  # no boundary for 'offsets'
  plower <- c(rep(-Inf,ncol(y)), pscale[(ncol(y)+1):(2*ncol(y))]/1e8) 
  if (is.null(pstart))
    pstart <- c(rep(0,ncol(y)), pscale[(ncol(y)+1):(2*ncol(y))])

  # if(verbose) 
  #   cat("pscale:", pscale, "\npstart:", pstart, "\nplower:", plower, "\n", sep="\t")
  
  control  <- list(trace=0, maxit=4000, parscale=pscale)
  optim.niter <- 10

  # a place to save the trajectory of estimated parameters along the iterations:
  params  <- matrix(NA, nrow=length(pscale), ncol=niter)

  sel <- rep(TRUE, nrow(y))
  for(lts.iter in 1:niter) {
    ## To save time, do in the beginning not fit with the full data set, but only
    ## with smaller random subsets. Later, go up to use full data.
    nrsel <- length(which(sel))
    nrsub <- max(1000, nrsel*lts.iter/niter) ## at least 1000
    if (nrsub < nrsel)
      sel[sel][runif(nrsel) > (nrsub/nrsel)] <- FALSE
    ## cat("lts.iter=", lts.iter, "nsel=", length(which(sel)), "\n")
    
    sy <- y[sel,]
    p0 <- pstart
    for (optim.iter in 1:optim.niter) {
      o  <- optim(par=p0, fn=llat, gr=grllat, method="L-BFGS-B",
                control=control, lower=plower)
      if (o$convergence==0) next
      
      cat("o$convergence=", o$convergence, ": ", o$message, "\n", sep="")
      if(o$convergence==52) {
        # ABNORMAL_TERMINATION_IN_LNSRCH
        # This seems to indicate that a stepwidth to go along the gradient could not be found,
        # probably because the start point p0 was already right at the optimum. Hence, try
        # again from a slightly different start point
        cat("lts.iter=", lts.iter, "optim.iter=", optim.iter, "pstart was", p0, "now trying ")
        p0 <- p0 + runif(length(pstart), min=0, max=0.01) * pscale 
        cat(p0, "\n")
      } else if(o$convergence==1) {
        # This seems to indicate that the max. number of iterations has been exceeded. Try again
        # with more
        cat("lts.iter=", lts.iter, "optim.iter=", optim.iter, "maxit was", control$maxit, "now trying ")
        control$maxit <- control$maxit*2
        cat(control$maxit, "\n")
      } else {
        stop(paste("Fatal error in optim, convergence=", o$convergence))  
      }
    }
    if (o$convergence!=0) 
      stop(paste("Fatal error in optim, not converged after", optim.niter, "times."))
    if (any(o$par[(ncol(sy)+1):(2*ncol(sy))]<0))
      stop(paste("Fatal error in optim, multiplicative normalization parameters have turned negative:", o$par))
    
    if(verbose)
      cat(sprintf("iter %2d: par=", as.integer(lts.iter)), sapply(o$par, function(x) sprintf("%9.3g",x)), "\n")
    
    # ----------------------------------------
    # selection of points in a LTS fashion
    # 1. calculate residuals; cf. llat()
    # ----------------------------------------
    offs <- matrix(o$par[    1  :  ncol(sy) ], nrow=nrow(y), ncol=ncol(sy), byrow=TRUE)
    facs <- matrix(o$par[(ncol(sy)+1):(2*ncol(sy))], nrow=nrow(y), ncol=ncol(sy), byrow=TRUE)
    asly <- asinh(offs + facs * y)         
    res  <- asly - rowMeans(asly)
    rsqres <- rowSums(res*res)
    hmean  <- rowSums(asly)
  
    # 2. select those data points within lts.quantile; do this separately for
    # each of nrslice slices along hmean
    nrslice  <- 5
    group    <- ceiling(rank(hmean)/length(hmean)*nrslice)
    grmed    <- tapply(rsqres, group, quantile, probs=lts.quantile)
    sel      <- rsqres <= grmed[group]
    
    params[,lts.iter] <- pstart <- o$par
  }
  hy = at(o$par)
  dimnames(hy) = dimnames(y)
  return(new("vsn.result", h=hy, params=params, sel=sel))
}

##-----------------------------------------------------------------
## .initvsn is called by .First.lib
##-----------------------------------------------------------------
.initvsn <- function(where) {
  ##------------------------------------------------------------
  ## define the class "vsn.result"
  ##------------------------------------------------------------
  setClass("vsn.result", where=where,
           representation(              ## its slots
                          h       = "matrix",
                          params  = "matrix",
                          sel     = "vector"), 
           prototype = list(            ## and default values
             h       = matrix(nrow=0, ncol=0),
             params  = matrix(nrow=0, ncol=0),
             sel     = logical(0)
             )
           )
  ##------------------------------------------------------------------
  ## params(vsn.result): will return the last column of matrix params
  ## i.e. the model parameters at the end of the iterations
  ##------------------------------------------------------------------
  if (!isGeneric("params"))
     setGeneric("params", function(object) standardGeneric("params"), where=where)
    
  setMethod("params", "vsn.result",
             function(object) object@params[,ncol(object@params)],
             where=where)
  
  ##------------------------------------------------------------
  ## plot method for vsn.result objects
  ##------------------------------------------------------------
  if (!isGeneric("plot"))
     setGeneric("plot", where=where, def=function(x, y, ...) standardGeneric("plot"))
  
  setMethod("plot", signature=c("vsn.result", "missing"), where=where,
     definition=function(x, ...) {
       if (ncol(x@h)<2 || ncol(x@h)*2 != nrow(x@params))
         stop("argument x is inconsistent")
       
       dots <- list(...)
       ind <- match("what", names(dots))
       if(!is.na(ind)) {
         what = dots$what
         dots <- dots[-ind]
       } else {
         what = "sdmean"
       }
       ind <- match("ranks", names(dots))
       if(!is.na(ind)) {
         ranks = dots$ranks
         dots <- dots[-ind]
       } else {
         ranks = TRUE
       }
       xlab <- ylab <- main <- ""
       ind <- match("xlab", names(dots))
       if(!is.na(ind)) {
         xlab <- dots$xlab
         dots <- dots[-ind]
       } 
       ind <- match("ylab", names(dots))
       if(!is.na(ind)) {
         ylab = dots$ylab
         dots <- dots[-ind]
       }
       ind <- match("main", names(dots))
       if(!is.na(ind)) {
         main = dots$main
         dots <- dots[-ind]
       }
       
       switch(what,
        sdmean = {
          cols <- c("black", "red")[as.numeric(x@sel)+1]
          n    <- length(x@sel)
          rkm  <- rank(rowMeans(x@h))
          sds  <- rowSds(x@h)
          ltsq <- length(which(x@sel))/n
          ## running quantile of width 2*dm
          dm        <- 0.1
          midpoints <- seq(dm, 1-dm, by=dm)
          rq.sds <- lapply(midpoints, function(mp) {
             median(sds[within(rkm/n, mp-dm, mp+dm)])
           })

          if(!is.logical(ranks))
            stop("argument ranks must be logical")
          if(ranks){
            px1 <- rkm; px2 <- midpoints*n
            if(xlab=="") xlab <- "rank of average intensity"
          } else {
            px1 <- rowMeans(x@h); px2 <- quantile(px1, probs=midpoints)
            if(xlab=="") xlab <- "average intensity"
          }
          if(ylab=="") ylab <- "standard deviation"
          
          args <- list(x=px1, y=sds, pch=".", col=cols, xlab=xlab, ylab=ylab)
          do.call("plot.default", append(args, dots))
          lines(px2, rq.sds, col="blue", type="b", pch=19)
        },
       ## 
       offsets = {
         d <- ncol(x@h)
         if(main=="") main <- what
         args <- list(x    = x@params[1,],
                      ylim = range(x@params[1:d,], na.rm=TRUE),
                      type = "b", pch=19, main=main, xlab=xlab, ylab=ylab)
         do.call("plot.default", append(args, dots)) 
         for (j in 2:d)
           lines(x@params[j,], type="b")
       },
       factors = {
         d <- ncol(x@h)
         if(main=="") main <- what
         args <- list(x    = x@params[d+1,],
                      ylim = range(x@params[d+(1:d),], na.rm=TRUE),
                      type = "b", pch=19, main=main, xlab=xlab, ylab=ylab)
         do.call("plot.default", append(args, dots)) 
         for (j in 2:d)
           lines(x@params[d+j,], type="b")
       },
       ##
       stop(paste("Unknown what=", what))
     ) ## end switch
   } ## end of function definition for "plot"
  )  ## end of setMethod("plot",...)
  
  ##------------------------------------------------------------
  ## print method for vsn.result objects
  ##------------------------------------------------------------
  if (!isGeneric("print"))
    setGeneric("print", where=where, def=function(x) standardGeneric("print"))
  
  setMethod("print", signature=c("vsn.result"), where=where,
     definition=function(x) {
       cat("vsn.result object\n",
         sprintf("%dx%d-matrix of transformed intensities\n", as.integer(nrow(x@h)), as.integer(ncol(x@h))),
         sprintf("%d-vector of transformation parameters\n", as.integer(nrow(x@params))))
     }
  )  ## end of setMethod("print")
  
  ##------------------------------------------------------------
  ## show method for vsn.result objects
  ##------------------------------------------------------------
  if (!isGeneric("show"))
    setGeneric("show", where=where, def=function(x) standardGeneric("show"))
  
  setMethod("show", signature=c("vsn.result"), where=where,
     definition=function(object) print(object) )
  
} ## end of .initvsn

##-----------------------------------------------------------------
## .First.lib: this function is called when the package is loaded
##-----------------------------------------------------------------
.First.lib <- function(lib, pkgname, where) {
  
  if(missing(where)) {
    where <- match(paste("package:", pkgname, sep=""), search())
    if(is.na(where)) {
      warning(paste("Not a package name: ",pkgname))
      return()
    }
    where <- pos.to.env(where)
  }
  .initvsn(where)
}


##------------------------------------------------------------
## a wrapper for vsn to be used as a normalization method in
## the package affy
##------------------------------------------------------------
normalize.Plob.vsn <- function(aplob, ...) {
  vsnres <- vsn(aplob@pm, niter=6)
  pm(aplob) = vsnh(aplob@pm, params(vsnres))
  mm(aplob) = vsnh(aplob@mm, params(vsnres))
  return(aplob)
}

## could also write "both" instead of "pm" if so desired
normalize.AffyBatch.vsn <- function(abatch, param) {
   if(!exists("indexProbes"))
     stop("Package affy must be loaded before calling normalize.AffyBatch.vsn")
   ind    <- unlist(indexProbes(abatch,"pm"))
   vsnres <- do.call("vsn", append(list(intensity(abatch)[ind,]), param))
   intensity(abatch) <- vsnh(intensity(abatch), params(vsnres))
   return(abatch)
}


##------------------------------------------------------------
## Some useful functions
## sqr   : square (ca. 10 times faster than ^2 !)
## rowSds: row standard deviations, cf. rowMeans  
## within: is x in interval [x1,x2] ?
##------------------------------------------------------------
   sqr <- function(x) { x*x }
rowSds <- function(x) { sqrt(rowSums(sqr(x-rowMeans(x)))/(ncol(x)-1)) }
within <- function(x, x1, x2) { x>=x1 & x<=x2 }

##---------------------------------------------------------
## Enumerate all the subsets of size k of the integers 1:n.
## The result is returned in a matrix with k rows and
## (n choose k) columns
## This function is not needed for VSN, but it is nice to 
## have it for the examples / the vignette and I haven't 
## yet found it anywhere else. 
##---------------------------------------------------------
nchoosek <- function(n, k) {
  if (!is.numeric(n)||!is.numeric(k)||is.na(n)||is.na(k)||length(n)!=1||length(k)!=1)
    stop("arguments must be non-NA numeric scalars.")
  if (k>n||k<0)
    stop("Arguments must satisfy 0 <= k <= n.")
  
  nck <- choose(n, k)
  res <- matrix(NA, nrow=k, ncol = nck)
  res[, 1] <- 1:k
  j <- 2
  repeat {
    res[, j] <- res[, j-1]
    i <- k
    repeat {
      res[i,j] <- res[i,j]+1
      if(res[i,j] <= n-(k-i))
        break
      i <- i-1
      stopifnot(i>=1)
    }
    if (i<k)
       res[(i+1):k,j] <- res[i,j] + 1:(k-i)
    j <- j+1
    if (j>nck) break
  }
  ## plausibility tests
  stopifnot(all(res[, nck]==(n-k+1):n))
  stopifnot(all(res<=n) && all(res>=1))
  return(res)
}

#-----------------------------------------------------------
# the "arsinh" transformation
#-----------------------------------------------------------
vsnh <- function(y, p) {
  if (!is.matrix(y) || !is.numeric(y))
    stop("vsnh: argument y must be a numeric matrix.\n")
  if (!is.vector(p) || !is.numeric(p) || any(is.na(p)))
    stop("vsnh: argument p must be a numeric vector with no NAs.\n")
  if (2*ncol(y) != length(p))
    stop("vsnh: argument p must be a vector of length 2*ncol(y).\n")

  offs <- matrix(p[         1  :  ncol(y) ], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  facs <- matrix(p[(ncol(y)+1):(2*ncol(y))], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  return(asinh(offs + facs * y))
}

#----------------------------------------------------------------------------------------
# debug: compare the analytical derivative with the numerical
#----------------------------------------------------------------------------------------
if (FALSE){
  cat('Checking  vsn.\n')
  if(!exists("vsnsample")) load("../data/vsnsample.RData") 

  tau <- vsnsample.cDNA
  z   <- vsn(tau, niter=0, verbose=TRUE)
  
  pscale <- c(rep(1, ncol(tau)), rep(1e-3, ncol(tau)))
  dp     <- pscale*0.0001
  nri    <- 10
  grs    <- array(NA, dim=c(2*ncol(tau), 2, nri))
  for (i in 1:nri) {
    p1   <- runif(4) * pscale
    amid      <- z$llat(p1)
    grs[,1,i] <- z$grllat(p1)
    for (j in 1:length(pscale)) {
      p2     <- p1
      p2[j]  <- p1[j]+dp[j]
      aright <- z$llat(p2)
      grs[j,2,i] <- (aright - amid)/dp[j]
    }
  }
  x11(width=6, height=6)
  par(mfrow=c(2,2))
  for (j in 1:4) {
    plot(grs[j,1,], grs[j,2,], cex=2, pch="x", main=paste(j),
         xlab="analytical", ylab="numerical")
    lines(1e10*c(-1,1), 1e10*c(-1,1), col="red")
  }
}

