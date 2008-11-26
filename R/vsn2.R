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

  o = .Call("vsn2_optim",
    v@x,
    p,
    istrat,
    v@reference@mu,
    v@reference@sigsq,
    v@optimpar,
    calibCharToInt(v@calib),
    PACKAGE="vsn")

  rv = new("vsn",
    coefficients = o$coefficients, 
    strata = v@strata,
    mu = o$mu,
    sigsq = o$sigsq,
    hx = matrix(NA_real_, nrow=nrow(v@x), ncol=0),
    lbfgsb = o$fail,
    hoffset = rep(NA_real_, nlevels(v@strata)),
    calib = v@calib)

  if (o$fail!=0L) {
    nrp = if(length(p)<=6L) length(p) else 6L
    msg = paste("** This is a diagnostic message on the performance of the optimizer,\n",
                "** it need not indicate a relevant problem:\n",
                sprintf("fail=%d\npstart[1:%d]=", o$fail, nrp),
                paste(signif(p[1:nrp], 4L), collapse=", "),
                sprintf("\n  coef[1:%d]=", nrp),
                paste(signif(coefficients(rv)[1:nrp], 4L), collapse=", "), "\n", sep="")
    ## warning(msg)
  }
  
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

  switch(vp@calib,
    affine = {
      if(length(vp@strata)>0L) {
        stopifnot(vp@ordered)
        istr = which(!duplicated(vp@strata))-1L
      } else {
        istr = 0L
      }
      stopifnot(length(istr)==nrs)

      istrat = integer(ncol(vp@x)*nrs+1L)
      for(i in 0:(ncol(vp@x)-1L))
        istrat[i*nrs + seq_len(nrs)] = i*nrow(vp@x) + istr

      istrat[length(istrat)] = nrow(vp@x)*ncol(vp@x)  ## end point
    },
         
    none = {
      if(nrs>1) stop("There must not be more than 1 stratum if 'calib' is 'none'")
      istrat = c(0L, nrow(vp@x)*ncol(vp@x))
    },
         
    stop(sprintf("Invalid value for 'vp@calib': %s", vp@calib))
  ) ## end of switch
  
  return(istrat)
}

progress = function(i, imax) {
  cat(if(i==0L) "\n" else "\r",
      sprintf("%3d", as.integer(i/imax*100)), "% done.", sep="")
}
 
##-------------------------------------------------------------------------
## LTS robust modification of the ML estimator
##-------------------------------------------------------------------------
vsnLTS = function(v) {

  ## for calculating a convergence criterion: earlier result
  oldhy   = Inf
  ## the number of iterations that have already met the convergence criterion
  cvgcCnt = 0
  
  ## integer version of "v@strata"
  intstrata = if(length(v@strata)==0L) rep(1L, nrow(v@x)) else as.integer(v@strata)
  facstrata = factor(intstrata, levels=paste(1:nlevels(v@strata)))
  stopifnot(!any(is.na(facstrata)))
    
  for(iter in seq_len(v@optimpar$cvg.niter)) {

    if(v@verbose)
      progress(iter-1, v@optimpar$cvg.niter)

    sv  = if(iter==1) v else v[whsel, ]
    rsv = vsnML(sv)

    ## if LTS.quantile is 1, then the following stuff is not necessary
    if (isSmall(v@lts.quantile-1))
      break
       
    ## apply to all data
    hy = vsn2trsf(x=v@x, p=coefficients(rsv), strata=intstrata, calib=v@calib)
    v@pstart = coefficients(rsv)

    ## Calculate residuals
    if(length(v@reference@mu)>0) {
      ## with reference:
      hmean = v@reference@mu
    } else {
      ## without reference:
      hmean = rowMeans(hy, na.rm=TRUE) 
      ## double check with what was computed in vsnML:
      if(iter==1L) {
        if(rsv@lbfgsb==0L) stopifnot(isSmall(rsv@mu-hmean))
      } else {
        if(rsv@lbfgsb==0L) stopifnot(isSmall(rsv@mu-hmean[whsel]))
        ## and create rsv@mu with mu of the right length (NA for the 'outliers' not in whsel)
        tmp = rep(NA_real_, nrow(v))
        tmp[whsel] = rsv@mu
        rsv@mu = tmp
      }
    }
    
    ## row variances
    rvar  = rowV(hy, mean=hmean, na.rm=TRUE)

    ## select those data points whose rvar is within the quantile; do this separately
    ## within each stratum, and also within 5 slices defined by hmean
    ## (see the SAGMB 2003 paper for details)
    nrslice = 5
    slice = ceiling(rank(hmean, na.last=TRUE)*(nrslice/length(hmean)))
    facslice = factor(slice, levels=paste(1:nrslice))
    stopifnot(!any(is.na(facslice)))
    
    grquantile = tapply(rvar, list(facslice, facstrata), quantile, probs=v@lts.quantile, na.rm=TRUE)
    if(any(is.na(grquantile)))
      stop(sprintf("Too many data points are NA (%d of %d), not enough data for fitting, am giving up.",
                   sum(is.na(sv@x)), length(sv@x)))

    ## Only use those datapoints with residuals less than grquantile,
    ##    but use all datapoints for slice 1 (the lowest intensity spots)
    whsel = which((rvar <= grquantile[ cbind(slice, intstrata) ]) |
                   slice == 1 )
    
    ## Convergence check
    ## after a suggestion from David Kreil, kreil@ebi.ac.uk
    if(v@optimpar$cvg.eps>0) {
      cvgc    = max(abs(hy - oldhy), na.rm=TRUE)
      cvgcCnt = if(cvgc<v@optimpar$cvg.eps) (cvgcCnt+1) else 0 
      if(cvgcCnt>=3)
        break
      oldhy = hy
    }

  } ## end of for-loop
  
  if(v@verbose){
    progress(1L, 1L)
    cat("\n")
  }
  
  return(rsv)
}

##----------------------------------------------------------------------------------
## vsnColumnByColumn
##----------------------------------------------------------------------------------
vsnColumnByColumn = function(v) {
  rlv = vsnLTS(v[,1L])
  d = dim(coefficients(rlv))
  n = ncol(v)
  stopifnot(length(d)==3, identical(d[2L], 1L))
  cf = array(NA_real_, dim=c(d[1L], n, d[3L]))
  cf[,1L,] = coefficients(rlv)
  for(j in seq_len(n)[-1L]) 
    cf[,j,] = coefficients(vsnLTS(v[,j]))

  return(new("vsn",
             coefficients = cf, 
             strata = v@strata,
             mu = v@reference@mu,
             sigsq = v@reference@sigsq,
             hoffset = rep(NA_real_,nlevels(v@strata)),
             calib = v@calib))
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

  res = if(length(v@reference@mu)>0L) {
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

  if( (v@subsample>0L) && (length(v@reference@mu)>0L) )
    stop("\n\nThe 'subsample' and 'normalization to reference' options cannot be mixed. ",
         "'normalization to reference' needs to use the same subsampling as the one ",
         "used for creating the reference. If you want to use 'normalization to reference', ",
         "please call this function without the 'subsample' argument.\n\n")

  wh = NULL
  
  if(v@subsample>0L) {
    if((length(v@strata)>0L) && (nlevels(v@strata)>1L)) { 
      sp = split(seq(along=v@strata), v@strata)
      wh = unlist(lapply(sp, sample, size=v@subsample))
    } else {
      wh = sample(nrow(v), size=v@subsample)
    }
  }

  if(length(v@reference@mu)>0L) {
    wh = which(!is.na(v@reference@mu))
  }
  

  if(!is.null(wh)){
    
    res = vsnStrata(v[wh, ])

    ## put back the results from subsampling 
    newmu = rep(NA_real_, nrow(v))
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
             calib        = "affine",
             pstart,
             minDataPointsPerStratum = 42L,
             optimpar   = list(),
             defaultpar = list(factr=5e7, pgtol=2e-4, maxit=60000L,
                               trace=0L, cvg.niter=7L, cvg.eps=0))
{
  storage.mode(x) = "double"
  storage.mode(subsample) = "integer"
  
  strata = if(missing(strata)){ 
    factor(integer(0), levels="all")
  } else {
    int2factor(strata)
  }
  
  stratasplit = if(length(strata)>0L) split(seq_len(nrow(x)), strata) else list(all=seq_len(nrow(x)))
  if(!(identical(names(stratasplit), levels(strata)) &&
       all(listLen(stratasplit)>minDataPointsPerStratum)))
    stop("One or more of the strata contain less than ", minDataPointsPerStratum, " elements.\n",
         "Please reduce the number of strata so that there is enough in each stratum.\n")
  
  if(missing(pstart))
    pstart = pstartHeuristic(x, stratasplit, calib)
  
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

  if(!(is.list(optimpar)&&all(names(optimpar)%in%optimparNames)&&!any(duplicated(names(optimpar)))))
    stop(paste("Names of elements of 'optimpar' must be from: ",
               paste("'", optimparNames, "'", collapse=", ", sep=""), ".", sep=""))
  opar = defaultpar
  opar[names(optimpar)] = optimpar
  
  v = new("vsnInput",
    x = x,
    strata = strata,
    pstart = pstart,
    reference = reference,
    lts.quantile = lts.quantile,
    optimpar = opar,
    subsample = subsample,
    verbose = verbose,
    calib = calib,
    ordered = FALSE)
  
  ## Print welcome message
  if (verbose)
    cat("vsn2: ", nrow(x), " x ", ncol(x), " matrix (", nlevels(strata), " strat",
        ifelse(nlevels(strata)==1L, "um", "a"), "). ", sep="")

  ## Here, the actual work is done.
  res = vsnSample(v)

  ## hoffset: compute from average scale factors between arrays, but separately for the different strata.
  ## coefficients 'cof' from are taken from the reference, if available.
  cof = coefficients( if(nrow(reference)==0L) res else reference )
  res@hoffset = log2(2*scalingFactorTransformation(rowMeans(cof[,,2L,drop=FALSE])))
  
  ## If necessary, calculate the data matrix transformed according to 'coefficients'
  if(returnData) {
    res@strata = strata
    res@hx = vsn2trsf(x = x,
      p = coefficients(res),
      strata = as.integer(res@strata),
      hoffset = res@hoffset,
      calib = calib)
  }
  
  stopifnot(validObject(res))

  if(verbose) {
    cat("Please use 'meanSdPlot' to verify the fit.\n")
  }
  return(res)
}

##---------------------------------------------------------------------
## The glog transformation
##--------------------------------------------------------------------
vsn2trsf = function(x, p, strata=numeric(0L), hoffset=NULL, calib="affine") {
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix.\n")
  
  if (!is.array(p) || !is.numeric(p)) stop("'p' must be a numeric array.\n")
  if (any(is.na(p))) stop("'p' must not contain NAs.\n")

  if(length(strata)==0L) {
    strata = rep(1L, nrow(x))
  } else {
    if(!is.integer(strata) || !is.vector(strata) || 
       length(strata)!=nrow(x) || any(is.na(strata)))
      stop("'strata' must be an integer vector of length nrow(x) with no NAs.")
  }
  nrstrata = max(strata)

  if(nrstrata==1L && length(dim(p))==2L)
    dim(p) = c(1L, dim(p))

  d2 = switch(calib,
    affine= ncol(x),
    none = 1,
    stop(sprintf("Invalid value of 'calib': %s", calib)))

  if(length(dim(p))!=3L ||
     dim(p)[1L] != nrstrata ||
     dim(p)[2L] != d2 ||
     dim(p)[3L] != 2L)
    stop("'p' has wrong dimensions.")
  
  hx = .Call("vsn2_trsf",
    x,
    as.vector(p),
    strata,
    calibCharToInt(calib),
    PACKAGE="vsn")
  
  dimnames(hx) = dimnames(x)
  
  ## see the 'value' section of the man page of 'vsn2'
  if(!is.null(hoffset)) {
    stopifnot(length(hoffset)==nrstrata)
    hx = hx/log(2) - hoffset[strata]
  }
  return(hx)
}

##-----------------------------------------------------------------------------
## The scaling factor transformation (function "f" in the vignette
## 'likelihoodcomputations.Rnw')
##------------------------------------------------------------------------------
scalingFactorTransformation = function(b) {
  .Call("vsn2_scalingFactorTransformation", b, PACKAGE="vsn")
}

##----------------------------------------------------------
## row-wise variances of a matrix (slightly more efficient
##  than geneplotter:rowVars since can pass 'mean')
##-----------------------------------------------------------
rowV = function(x, mean, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<1]  = NA
  if(missing(mean))
    mean=rowMeans(x, ...)
  return(rowSums(sqr(x-mean), ...)/(n-1))
}

##--------------------------------------------------------------
## A heuristic to set the start parameters for offset and scale.
##--------------------------------------------------------------
pstartHeuristic = function(x, sp, calib) {
  
  d2 = switch(calib,
    affine= ncol(x),
    none = 1,
    stop(sprintf("Invalid value of 'calib': %s", calib)))

  pstart = array(0, dim=c(length(sp), d2, 2))
  pstart[,,2] = 1
  return(pstart)
}

##-------------------- 
## integer to factor
##-------------------- 
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

##-------------------------------------------------------
## check if all elements of a vector are close to 0
##--------------------------------------------------------
isSmall = function(x, tol=sqrt(.Machine$double.eps)) (max(abs(x))<tol)

