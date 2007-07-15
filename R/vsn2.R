##----------------------------------------------------------------------
## Robust calibration and variance stabilization
##   revised version vsn2
## (C) Wolfgang Huber <huber@ebi.ac.uk> 2002-2007
## contributions from Markus Ruschhaupt, Dennis Kostka, David Kreil
##
## The computations are organised as a set of nested function calls.
## The innermost is vsnML, which calls the C code.
## The next one is vsnLTS, which does the robust modification (least
##   trimmed sum of squares, LTS) of the maximum likelihood (ML) estimator.
## The next one is vsnStrata, which if necessary reorders the data so
##   that features from the same stratum are adjacent.
## The next one is vsnSample, which allows for sampling.
## The next one is vsnMatrix, a function that accepts a matrix with the
##   data and the various user-defined parameters
## The outmost ones are various S4 methods that call vsnMatrix.
##----------------------------------------------------------------------

##----------------------------------------------------------------------------
## ML estimator
##-----------------------------------------------------------------------------
vsnML = function(v) {

  p = as.vector(v@pstart)
  istrat = calcistrat(v) ## pointers to the starts of strata

  o = .Call("vsn2_optim", v@x, p, istrat, v@reference@mu,
    v@reference@sigsq, v@optimpar, PACKAGE="vsn")

  rv = new("vsn", coefficients=o$coefficients, 
           mu=o$mu, sigsq=o$sigsq,
           strata=v@strata)
  
  if (o$fail!=0L) {
    nrp = if(length(p)<=6L) length(p) else 6L
    msg = paste(sprintf("fail=%d\npstart[1:%d]=", o$fail, nrp),
                paste(signif(p[1:nrp], 4L), collapse=", "),
                sprintf("\n  coef[1:%d]=", nrp),
                paste(signif(coefficients(rv)[1:nrp], 4L), collapse=", "), sep="")
    warning(msg)
  }
  
  if (any(coefficients(rv)[,,2L]<0))
    stop("Internal error: likelihood optimization produced negative scale factor estimates.\n",
         "Please contact the package maintainer.\n")
  
  return(rv)
}

##------------------------------------------------------------
## istrat[j] is the starting position for stratum j
## i.e. istrat[j]...istrat[j+1] are the elements of ysel that
## belong to stratum j (using C indexing convention,
## i.e. starting at 0).
## Note: the counting over the different samples is folded into j,
## i.e., if there are 4 strata on the array and 3 samples, then
## j runs from 1:12
##------------------------------------------------------------
calcistrat = function(vp) {
  nrs = nlevels(vp@strata)
  one = as.integer(1)

  if(length(vp@strata)>0) {
    stopifnot(vp@ordered)
    istr = which(!duplicated(vp@strata))-one
  } else {
    istr = as.integer(0)
  }
  stopifnot(length(istr)==nrs)

  istrat = integer(ncol(vp@x)*nrs+one)
  for(i in 0:(ncol(vp@x)-one))
    istrat[i*nrs + seq_len(nrs)] = i*nrow(vp@x) + istr

  istrat[length(istrat)] = nrow(vp@x)*ncol(vp@x)  ## end point
  return(istrat)
}

progress = function(i, imax) {
  if(i>0)
    cat("\b\b\b\b\b\b\b\b\b\b")
  cat(sprintf("%3d", as.integer((1-((i-imax)/imax)^2)*100)), "% done.", sep="") ## 10 characters
}
 
##-------------------------------------------------------------------------
## LTS robust modification of the ML estimator
##-------------------------------------------------------------------------
vsnLTS = function(v) {

  ## for calculating a convergence criterion: earlier result
  oldhy   = Inf
  
  ## the number of iterations that have already met the convergence criterion
  cvgcCnt = 0L

  ## integer version of "v@strata"
  intStrat = if(length(v@strata)==0L) rep(1L, nrow(v@x)) else as.integer(v@strata)

  if(v@verbose)
      progress(0, v@optimpar$cvg.niter)

  for(iter in seq_len(v@optimpar$cvg.niter)) {
    sv  = if(iter==1) v else v[whsel, ]
    rsv = vsnML(sv)

    ## apply to all data
    hy = vsn2trsf(v@x, coefficients(rsv), intStrat)
    v@pstart = coefficients(rsv)

    ## if LTS.quantile is 1, then the following stuff is not necessary
    if (abs(v@lts.quantile-1)<sqrt(.Machine$double.eps))
      break
       
    ## Calculate residuals
    hmean  = if(length(v@reference@mu)>0) {
      v@reference@mu           ## WITH reference
    } else {
      rowMeans(hy, na.rm=TRUE)   ## WITHOUT reference
    }

    ## row variances
    rsd  = rowVars(hy, mean=hmean, na.rm=TRUE)

    ## select those data points whose rsd is within the quantile; do this separately
    ## within each stratum, and also within strata defined by hmean
    ## (see the SAGMB 2003 paper for details)
    nrslice = 5
    slice   = ceiling(rank(hmean, na.last=TRUE)/length(hmean)*nrslice)
    slice   = factor((intStrat-1)*nrslice + slice)
    grmed   = tapply(rsd, slice, quantile, probs=v@lts.quantile, na.rm=TRUE)
    if(any(is.na(grmed)))
      stop(sprintf("Too many data points are NA (%d of %d), not enough data for fitting, am giving up.",
                   sum(is.na(sv@x)), length(sv@x)))
    
    meds    = grmed[as.character(slice)]
    whsel   = which(rsd <= meds)

    ## diagnostic plot (most useful for d=2)
    ## plot(hy, pch=".", col=2-sel, main=sprintf("iter=%d", iter))
    
    ## Convergence check
    ## after a suggestion from David Kreil, kreil@ebi.ac.uk
    if(v@optimpar$cvg.eps>0) {
      cvgc    = max(abs(hy - oldhy), na.rm=TRUE)
      cvgcCnt = if(cvgc<v@optimpar$cvg.eps) (cvgcCnt+1) else 0 
      if(cvgcCnt>=3)
        break
      oldhy = hy
    }

    if(v@verbose)
      progress(iter, v@optimpar$cvg.niter)

  } ## end of for-loop
  
  if(v@verbose) {
    progress(v@optimpar$cvg.niter, v@optimpar$cvg.niter)
    cat("\n")
  }
  return(rsv)
}

##----------------------------------------------------------------------------------
## vsnColumnByColumn
##----------------------------------------------------------------------------------
vsnColumnByColumn = function(v) {
  rlv = vsnLTS(v[,1])
  d = dim(coefficients(rlv))
  n = ncol(v)
  stopifnot(d[2]==1L)
  cf = array(as.numeric(NA), dim=c(d[1], n, d[3]))
  cf[,1,] = coefficients(rlv)
  if(n>1) {
    for(j in 2:n) {
      cf[,j,] = coefficients(vsnLTS(v[,j]))
    }
  } 
  return(new("vsn", coefficients=cf, 
             mu = v@reference@mu, sigsq = v@reference@sigsq,
             strata=v@strata))
}

##----------------------------------------------------------------------------------
## vsnStrata: if necessary,
## reorder the rows of the matrix so that each stratum sits in a contiguous block
##----------------------------------------------------------------------------------
vsnStrata = function(v) {

  ord = NULL
  if(nlevels(v@strata)>1) {
    ord = order(v@strata)
    ## reorder the rows of the matrix so that each stratum sits in a contiguous block
    v = v[ord,]
  }
  v@ordered = TRUE

  res = if(length(v@reference@mu)>0) {
    vsnColumnByColumn(v)
  } else {
    vsnLTS(v)
  }

  ## reverse the ordering
  res@mu[ord] = res@mu
  
  return(res)
} 

##----------------------------------------------------------------------------------
## vsnSample: allows for (balanced) subsetting
##----------------------------------------------------------------------------------
vsnSample = function(v) {

  if(v@subsample>0) {
    if((length(v@strata)>0) && (nlevels(v@strata)>1)) { 
      sp = split(seq(along=v@strata), v@strata)
      wh = unlist(lapply(sp, sample, size=v@subsample))
    } else {
      wh = sample(nrow(v), size=v@subsample)
    }
    res = vsnStrata(v[wh, ])

    ## put back the results from subsampling 
    newmu = numeric(nrow(v))
    newmu[wh] = res@mu
    res@mu = newmu
    
  } else {
    res = vsnStrata(v)
  }  
  return(res)
}

##----------------------------------------------------------------------
## vsnMatrix
##----------------------------------------------------------------------
vsnMatrix = function(x,
  reference,
  strata,
  lts.quantile = 0.9,
  subsample    = 0L,
  verbose      = interactive(),
  returnData   = TRUE,
  pstart,
  optimpar   = list(),
  defaultpar = list(factr=5e7, pgtol=0, lower=2e-5, maxit=60000L, trace=0L, cvg.niter=1L, cvg.eps=0)) {

  storage.mode(x) = "double"
  storage.mode(subsample) = "integer"
  
  strata = if(missing(strata)){ 
    factor(integer(0), levels="all")
  } else {
    int2factor(strata)
  }
  
  minDataPointsPerStratum = 42L
  stratasplit = if(length(strata)>0L) split(seq_len(nrow(x)), strata) else list(all=seq_len(nrow(x)))
  if(!(identical(names(stratasplit), levels(strata)) &&
       all(listLen(stratasplit)>minDataPointsPerStratum)))
    stop("One or more of the strata contain less than ", minDataPointsPerStratum, " elements.\n",
         "Please reduce the number of strata so that there is enough in each stratum.\n")
  
  if(missing(pstart))
    pstart = pstartHeuristic(x, stratasplit)
  
  if(missing(reference)) {
    if(ncol(x)<=1L)
      stop("'x' needs to have 2 or more columns if no 'reference' is specified.")
    reference = new("vsn")
  } else {
    if(nrow(reference)!=nrow(x))
      stop("'nrow(reference)' must be equal to 'nrow(x)'.")
    if(nrow(reference)!=length(reference@mu))
      stop(sprintf("The slot 'reference@mu' has length %d, but expected is n=%d", nrow(reference)))
  }

  if(!(is.list(optimpar)&&all(names(optimpar)%in%names(defaultpar))))
    stop("'optimpar' must be a list whose elements have the same names as elements of 'defaultpar'.")
  opar = defaultpar
  opar[names(optimpar)] = optimpar
  
  v = new("vsnInput",
    x      = x,
    strata = strata,
    pstart = pstart,
    reference = reference,
    lts.quantile = lts.quantile,
    optimpar = opar,
    subsample = subsample,
    verbose   = verbose,
    ordered   = FALSE)
  
  ## Print welcome message
  if (verbose)
    cat("vsn: ", nrow(x), " x ", ncol(x), " matrix (", nlevels(strata), " strat",
        ifelse(nlevels(strata)==1L, "um", "a"), "). ", sep="")

  res = vsnSample(v)

  ## If necessary, calculate the data matrix transformed according to 'coefficients'
  ## Apply an irrelevant affine transformation to make users happy (with coefficients
  ## from reference, if applicable).
  if(returnData) {
    res@strata=strata
    trsfx = vsn2trsf(x, coefficients(res), strata=
      if(length(strata)==0L) rep(1L, nrow(x)) else as.integer(strata))
    res@hx = trsf2log2scale(trsfx,  coefficients=
      if(nrow(reference)==0L) coefficients(res) else coefficients(reference))
  }
  
  return(res)
}

##---------------------------------------------------------------------
## trsf2log2scale
##--------------------------------------------------------------------
trsf2log2scale = function(x, coefficients)
  (x-mean(log(2/coefficients[,,2])))/log(2)

##---------------------------------------------------------------------
## The glog transformation
##--------------------------------------------------------------------
vsn2trsf = function(x, p, strata) {
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix.\n")
  
  if (!is.array(p) || !is.numeric(p) || any(is.na(p)))
    stop("'p' must be an array with no NAs.\n")

  if(missing(strata)) {
    strata = rep(1L, nrow(x))
  } else {
    if(!is.integer(strata) || !is.vector(strata) || 
       length(strata)!=nrow(x) || any(is.na(strata)))
      stop("'strata' must be an integer vector of length nrow(x) with no NAs.")
  }
  nrstrata = max(strata)
  
  if(nrstrata==1L && length(dim(p))==2L)
    dim(p) = c(1L, dim(p))
  
  if(length(dim(p))!=3L || dim(p)[1L]!=nrstrata || dim(p)[2L]!=ncol(x) || dim(p)[3L]!=2L)
    stop("'p' has wrong dimensions.")
  if (any(p[,,2L]<=0))
    stop("'p' contains invalid values: factors must be non-negative.")

  hx = .Call("vsn2_trsf", x, as.vector(p), strata, PACKAGE="vsn")

  dimnames(hx) = dimnames(x)
  return(hx)
}

##----------------------------------------------------------
## helper function: row-wise variances of a matrix
##-----------------------------------------------------------

rowVars = function(x, mean, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<1] = NA
  if(missing(mean))
    mean=rowMeans(x, ...)
  return(rowSums(sqr(x-mean), ...)/n)
}

##
## A heuristic to set the start parameters for offset and scale.
pstartHeuristic = function(x, sp) {
  pstart = array(0, dim=c(length(sp), ncol(x), 2L))
#  for(i in seq_along(sp)) {
#    for(j in seq_len(ncol(x))) {
#      rg = quantile(x[sp[[i]], j], probs=c(0.10, 0.75), na.rm=TRUE)
#      z = 8/(rg[2]-rg[1])
#      pstart[i,j,] = c(2-rg[1]*z, z)
#    }
#  }
  pstart[,,2] = 1
  return(pstart)
}

## -------------------- 
## integer to factor
## -------------------- 
int2factor = function(strata) {
  if(is.factor(strata))
    return(strata)
  if(!is.integer(strata))
    stop("'strata' needs to be an integer vector.")
  ssu = sort(unique(strata))
  if(!identical(ssu, seq(along=ssu)))
    stop("'strata' needs to be an integer vector whose values cover ",
         "the range of  1...n.")
  factor(strata, levels=ssu)
}
