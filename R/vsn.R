##-----------------------------------------------------------------
## Robust calibration and variance stabilization
## (C) Wolfgang Huber 2002-2003
## w.huber@dkfz.de
##-----------------------------------------------------------------
require(Biobase) || stop("Cannot load without package \"Biobase\"")

##------------------------------------------------------------
## vsn: the main function of this library
##------------------------------------------------------------
vsn = function(intensities, lts.quantile=0.5, niter=10, verbose=TRUE, pstart=NULL) {
  ## Bureaucracy step 1: make sure are the arguments are valid and plausible
  mess =  badarg = NULL
  if (!is.numeric(lts.quantile) || (length(lts.quantile)!=1) ||
     (lts.quantile<0.5) || (lts.quantile>1)) {
    badarg = "lts.quantile"
    mess   = "Please specify a scalar between 0.5 and 1."
  }
  if (!is.numeric(niter) || (length(niter)!=1) || (niter<1)) {
    badarg = niter
    mess   = "Please specify a number >=1"
  }
  if (!is.null(pstart))
    if (!is.numeric(pstart) || length(pstart) != 2*d || any(is.na(pstart)) ||
        any(pstart[(ncol(intensities)+1):(2*ncol(intensities))]<=0)) {
      badarg = "pstart"
      mess   = paste("Please specify a numeric vector of length", 2*ncol(intensities))
    }
  if (!is.logical(verbose)) {
    badarg = "verbose"
    mess   = "Please specify a logical value"
  }

  ## Bureaucracy step 2: extract the intensity matrix from the argument "intensities"
  y = switch(class(intensities),
     matrix     = {  if (!is.numeric(intensities)) {
                       badarg = "intensities"
                       mess   = "Please specify a numeric matrix"
                     }
                     intensities
                   },
     data.frame = {  if (!all(sapply(intensities, is.numeric))) {
                       badarg = "intensities"
                       mess   = "Please specify a data.frame that has only numeric columns"
                     }
                     as.matrix(intensities)
                  },
     exprSet    = { tmp <- exprs(intensities)
                    colnames(tmp) <- rownames(pData(intensities))
                    tmp
                  },
     marrayRaw  = { nrslides <- ncol(intensities@maRf)
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
    { badarg = "intensities"
      mess   = paste("It has class ", class(intensities),
                      ". Permitted are: matrix, data.frame, exprSet, marrayRaw", sep="")
    }
    )  ## end of switch statement

  if (any(is.na(y))) {
    badarg = "intensities"
    mess   = paste("It must not contain NA values.\n",
             "This could indicate that the input data has already undergone some\n",
             "thresholding or transformation (log?), and may not satisfy the\n",
             "requirements of the multiplicative-additive noise model.\n",
             "Otherwise, consider calling vsn on a subset of data where all values\n",
             "are defined, and then use vsnh.\n")
  }

  ## Error handling
  if(!is.null(mess)) {
    if(badarg=="intensities") {
      mess = paste("The argument ", badarg, " has an invalid value.\n", mess, "\n", sep="")
    } else {
      mess = paste("The argument ", badarg, " has an invalid value:", get(badarg), "\n", mess, "\n", sep="")
    }
    stop(mess)
  }

  ## Print welcome message
  if (verbose)
    cat("vsn is working on ", nrow(y), " x ", ncol(y), " matrix, with lts.quantile=", signif(lts.quantile, 2),
        "; please wait for ", niter+1, " dots:\n.", sep="")

  ##----------------------------------------------------------
  ## Profile log likelihood of the model
  ##
  ##    asinh(a_i + b_i * y_ki) = m_k + epsilon_ki
  ##
  ## where k=1..n, i=1..d, y is a n-by-d matrix,
  ## (a_i, b_i) are real parameters to be estimated
  ## epsilon_ki are i.i.d N(0, sigma^2)
  ## and sigma, and all m_k are estimated from the data
  ## ("profiling")
  ##
  ## Argments:
  ## p numeric vector of length 2*d, with elements a_1,...,a_d, b_1,...,b_d
  ##
  ## Variables stored in environment "ws":
  ## sy     the matrix of the y_ki
  ## ly     the matrix of the a_i + b_i * y_ki
  ## asly   the matrix of the asinh(a_i + b_i * y_ki)
  ## res    the matrix of the residuals epsilon_ki
  ## nrs    number of rows of the above matrices
  ## ncs    number of columns of the above matrices
  ## ssq    sum of squared residuals sum(res^2)
  ## p      a place where ll() stores its call arguments. From this, when grll()
  ##        is called it can double-check whether it is indeed also called
  ##        with the same arguments (see below for details)
  ##
  ## Return value:
  ## Profile Log-Likelihood
  ##----------------------------------------------------------
  ll = function(p) {
    assign("p", p, envir=ws)
    with(ws, {
      offs = p[ 1:ncs     ]
      facs = p[(1:ncs)+ncs]

      ## sy is an nxd matrix, offs and facs are d-vectors
      ## the cycling goes column-wise, hence transpose!
      ly   = offs + facs * t(sy)
      asly = t(asinh(ly))

      ## residuals
      res = asly - rowMeans(asly)
      ssq = sum(res*res)
      rv  = nrs*ncs/2*log(ssq) - sum(log(facs/sqrt(1+ly*ly)))
    } )
    return(get("rv", ws))
  }
  ##--------------------------------------------------------------
  ## Gradient of the profile log likelihood
  ##--------------------------------------------------------------
  grll <- function(p) {
    ## Generally, optim() will call the gradient of the objective function (gr)
    ## immediately after a call to the objective (fn) at the same parameter values.
    ## Anyway, we like to doublecheck
    if(any(p!=get("p", ws))) {
      mess = paste(
        "\n\n\The function grll (likelihood-gradient) was called with different\n",
        "parameters than the previous call of ll (likelihood-value).\n",
        "This should never happen. Please contact package maintainer at\n",
        "w.huber@dkfz.de\n\n")
      error(mess)
    }
    with(ws, {
      dhda      = 1/sqrt(1+ly*ly)
      dlndhdyda = -ly/(1+ly*ly)
      gra       = nrs*ncs/ssq*res*t(dhda) - t(dlndhdyda)
      rv        = c(colSums(gra), colSums(gra*sy) - nrs/p[(ncs+1):(ncs*2)])
    } )
    return(get("rv", ws))
  }
  ##----------------------------------------------------------------------
  ## guess a parameter scale, set boundaries for optimization,
  ## and, if they are not user-supplied, set the start parameters
  ##----------------------------------------------------------------------
  pscale  <- plower <- numeric(2*ncol(y))
  pscale[1:ncol(y)] <- 1
  for (j in 1:ncol(y)) {
    pscale[ncol(y)+j] <- 1/diff(quantile(y[,j], probs=c(0.25, 0.75)))
  }
  ## lower boundary for 'factors': a small positive value
  ## no boundary for 'offsets'
  plower <- c(rep(-Inf,ncol(y)), pscale[(ncol(y)+1):(2*ncol(y))]/1e8)
  if (is.null(pstart))
    pstart <- c(rep(0,ncol(y)), pscale[(ncol(y)+1):(2*ncol(y))])

  ## factr controls the convergence of the "L-BFGS-B" method. Convergence
  ## occurs when the reduction in the objective is within this factor of
  ## the machine tolerance. Default is 1e7, that is a tolerance of about
  ## 1e-8. Here we use 5e8 to save a little time.
  control     <- list(trace=0, maxit=4000, parscale=pscale, factr=5e8)
  optim.niter <- 10

  ## a place to save the trajectory of estimated parameters along the iterations:
  params  <- matrix(NA, nrow=length(pscale), ncol=niter)

  ## Workspace. This has two purposes: 1. we don't need to pass y or sy around,
  ## 2. the intermediate results of "ll" can be used by "grll"
  ws = new.env(hash=TRUE)
  assign("timell",   numeric(0), envir=ws)
  assign("timegrll", numeric(0), envir=ws)

  sel <- rep(TRUE, nrow(y))
  for(lts.iter in 1:niter) {
    assign("sy",  y[sel,],            envir=ws)
    assign("nrs", length(which(sel)), envir=ws)
    assign("ncs", ncol(y),            envir=ws)

    p0 <- pstart
    for (optim.iter in 1:optim.niter) {
      o  <- optim(par=p0, fn=ll, gr=grll, method="L-BFGS-B",
                  control=control, lower=plower)
      if (o$convergence==0) next

      if(o$convergence==52) {
        ## ABNORMAL_TERMINATION_IN_LNSRCH
        ## This seems to indicate that a stepwidth to go along the gradient could not be found,
        ## probably because the start point p0 was already right at the optimum. Hence, try
        ## again from a slightly different start point
        ## cat("lts.iter=", lts.iter, "optim.iter=", optim.iter, "pstart was", p0, "now trying ")
        p0 <- p0 + runif(length(pstart), min=0, max=0.01) * pscale
        ## cat(p0, "\n")
      } else if(o$convergence==1) {
        ## This seems to indicate that the max. number of iterations has been exceeded. Try again
        ## with more
        ## cat("lts.iter=", lts.iter, "optim.iter=", optim.iter, "maxit was", control$maxit, "now trying ")
        control$maxit <- control$maxit*2
        ## cat(control$maxit, "\n")
      } else {
        stop(paste("Likelihood optimization: the function optim() returned the value convergence=",
                   o$convergence, "\nPlease make sure your data is good.",
                   "If so, contact the package maintainer.\n", sep=""))
      }
    }

    mess = NULL
    if (o$convergence!=0)
      mess = paste("Likelihood optimization did not converge even after", optim.niter, "calls to optim().",
                   "\nPlease make sure your data is good. If the problem persists,",
                   "\nplease contact the package maintainer.\n")

    if (any(o$par[(ncol(y)+1):(2*ncol(y))]<0))
      mess = paste("Likelihood optimization produced negative parameter estimates in spite of constraints.",
                   "\nPlease contact the package maintainer.\n")

    if(!is.null(mess))
      stop(mess)

    if(verbose)
      cat(".")

    # ----------------------------------------
    # selection of points in a LTS fashion
    # 1. calculate residuals; cf. ll()
    # ----------------------------------------
    offs <- matrix(o$par[    1      :   ncol(y) ], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
    facs <- matrix(o$par[(ncol(y)+1):(2*ncol(y))], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
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
  if(verbose)
    cat("\n")

  ## Profiling results:
  ##  cat3 = function(x)
  ##    paste(length(x), " calls, min time=", min(x), "  max time=", max(x),
  ##          "  mean=", signif(mean(x),3), "\n", sep="")
  ##  cat("  ll:", cat3(get("timell",   ws)))
  ##  cat("grll:", cat3(get("timegrll", ws)))

  return(new("vsn.result", h=vsnh(y, o$par), params=params, sel=sel))
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

  setMethod("show", signature=c("vsn.result"), where=where,
     definition=function(object) print(object) )

} ## end of .initvsn


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

##---------------------------------------------------------------------
## the "arsinh" transformation
## Note: an overall additive constant has no effect on the Delta-h, i.e.
## on the differences between transformed intensities. We choose the
## additive constant such that h_1(y) \approx log(y) for y\to\infty.
##---------------------------------------------------------------------
vsnh <- function(y, p) {
  if (!is.matrix(y) || !is.numeric(y))
    stop("vsnh: argument y must be a numeric matrix.\n")
  if (!is.vector(p) || !is.numeric(p) || any(is.na(p)))
    stop("vsnh: argument p must be a numeric vector with no NAs.\n")
  if (2*ncol(y) != length(p))
    stop("vsnh: argument p must be a vector of length 2*ncol(y).\n")

  offs <- matrix(p[         1  :  ncol(y) ], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  facs <- matrix(p[(ncol(y)+1):(2*ncol(y))], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  hy   <- asinh(offs + facs * y) - log(2*facs[1])
  dimnames(hy) = dimnames(y)
  return(hy)
}
