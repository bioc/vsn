##----------------------------------------------------------------------
## Robust calibration and variance stabilization
## (C) Wolfgang Huber <huber@ebi.ac.uk> 2002-2004
## With contributions from Markus Ruschhaupt, Dennis Kostka, David Kreil
##----------------------------------------------------------------------
##----------------------------------------------------------------------
## vsn: the main function of this library
##----------------------------------------------------------------------
vsn <- function(intensities,
                lts.quantile = 0.5,
                verbose      = TRUE,
                niter        = 10,
                cvg.check    = NULL,
                describe.preprocessing=TRUE,
                pstart,
                strata) {
  y <- getIntensityMatrix(intensities, verbose)
  d <- ncol(y)

  ## Make sure the arguments are valid and plausible
  if (!is.numeric(lts.quantile) || (length(lts.quantile)!=1) ||
    (lts.quantile<0.5) || (lts.quantile>1))
      stop("'lts.quantile' must be a scalar between 0.5 and 1.")
  if (!is.numeric(niter) || (length(niter)!=1) || (niter<1))
    stop("'niter' must be a number >=1.")
  if (!is.logical(verbose))
    stop("'verbose' must be a logical value.")

  if(missing(strata)) {
    strata <- rep(as.integer(1), nrow(y))
  } else {
    if(!is.integer(strata) || !is.vector(strata) ||
       length(strata)!=nrow(y) || any(is.na(strata)))
      stop("'strata' must be an integer vector of length nrow(y) with no NAs.")
  }

  nrstrata <- max(strata)
  if(missing(pstart)) {
    pstart <- array(0, dim=c(nrstrata, d, 2))
    for(i in 1:d)
      pstart[,i,2] <- 1/diff(quantile(y[,i], probs=c(0.25, 0.75)))
  } else {
    if(!is.array(pstart) || length(dim(pstart))!=3)
      stop("'pstart' must be a 3D array.")
    if(!all(dim(pstart)==c(nrstrata, d, 2)))
      stop(paste("dimensions of 'pstart' do not match. They should be ",
        paste(nrstrata, d, 2, sep=" x "), ", but are ",
        paste(dim(pstart), collapse=" x "), ".", sep=""))
  }

  minperstratum <- as.integer(42/lts.quantile)
  sstr <- sum(table(strata) < minperstratum)
  if(sstr>0) {
    mess <- paste("*** There are less than", minperstratum, "data points in", sstr,
      "of the strata.\n*** The fitted parameters may be unreliable.\n")
    if(lts.quantile<0.9)
      mess <- paste(mess, "*** You could try to increase the value of 'lts.quantile'.\n", sep="")
    if(nrstrata>1)
      mess <- paste(mess, "*** You could try to reduce the number of strata.\n", sep="")
    warning(mess)
  }

  ## reorder the rows of the matrix so that each stratum sits in a contiguous block
  ordstrata <- order(strata)
  reord     <- order(ordstrata)
  y         <- y[ordstrata,]
  strata    <- strata[ordstrata]
  optim.niter <- 10

  ## Print welcome message
  if (verbose)
    cat("vsn: ", nrow(y), " x ", d, " matrix (", nrstrata, " strat",
        ifelse(nrstrata==1, "um", "a"), "). Please wait for ",
        niter, " dots: ", sep="")

  ## a place to save the trajectory of estimated parameters along the iterations:
  params <- array(NA, dim=c(dim(pstart), niter))

  ##--------------------------------------------------
  ## begin of the outer LL iteration loop
  ##--------------------------------------------------
  oldhy   <- Inf  ## for calculating a convergence criterion: earlier result
  cvgcCnt <- 0    ## counts the number of iterations that have already met the convergence criterion
  sel     <- rep(TRUE, nrow(y))
  for(lts.iter in 1:niter) {
    ysel   <- y[sel,]
    strsel <- strata[sel]

    ## istrat[j] is the starting position for stratum j
    ## i.e. istrat[j]...istrat[j+1] are the elements of ysel that
    ## belong to stratum j (using C indexing convention,
    ## i.e. starting at 0).
    ## Note: the counting over the different samples is folded into j,
    ## i.e., if there are 4 strata on the array and 3 samples, then
    ## j runs from 1:12
    istr <- which(!duplicated(strsel))-1
    stopifnot(length(istr)==nrstrata, all(diff(strsel[istr])>0))

    istrat <- numeric(d*nrstrata+1)
    istrat[d*nrstrata+1] <- length(ysel)  ## end point
    for(i in 1:d)
      istrat[(i-1)*nrstrata + (1:nrstrata)] <- ((i-1)*nrow(ysel) + istr)
    istrat <- as.integer(istrat)

    p0 = pstart
    for (optim.iter in 1:optim.niter) {
      ## cat(lts.iter, ":", optim.iter, "\t", paste(signif(p0, 4), collapse=" "), "\n", sep="")
      optres <- .Call("vsnc", ysel, as.vector(p0), istrat, TRUE, PACKAGE="vsn")
      stopifnot(length(optres)==2*nrstrata*d+1)
      conv <- round(optres[length(optres)])
      par  <- array(optres[-length(optres)], dim=dim(p0))
      if (conv==0)
        break
      if(conv==52) {
        ## ABNORMAL_TERMINATION_IN_LNSRCH
        ## This seems to indicate that a stepwidth to go along the gradient could not be found,
        ## probably because the start point p0 was already too close to the optimum (?). Hence, try
        ## again from a slightly different start point
        p0 = p0 + runif(length(pstart), min=0, max=0.01)
        if(verbose)
          cat("(CONV=52, restarting with new initial parameters)")
      } else {
        stop(paste("Likelihood optimization: the function optim() returned the value convergence=",
                   conv, "\nPlease make sure your data is good. ",
                   "If so, contact the package maintainer.\n", sep=""))
      }
    }

    if (conv!=0)
      stop(paste("Likelihood optimization did not converge even after", optim.niter, "calls to optim().",
           "\nPlease make sure your data is good. If the problem persists,",
           "\nplease contact the package maintainer.\n"))
    if (any(par[,,2]<0))
      stop(paste("Likelihood optimization produced negative parameter estimates in spite of constraints.",
           "\nPlease contact the package maintainer.\n"))

    if(verbose)
      cat(".")

    ## selection of points in a LTS fashion:
    ## calculate residuals
    hy     <- vsnh(y, par, strata)
    hmean  <- rowMeans(hy)
    sqres  <- hy - hmean
    sqres  <- rowSums(sqres*sqres) ## squared residuals

    ## select those data points within lts.quantile; do this separately
    ## within each stratum, and also within strata defined by hmean
    ## (see the SAGMB 2003 paper for details)
    nrslice <- 5
    group   <- ceiling(rank(hmean)/length(hmean)*nrslice)
    group   <- factor((strata-1)*nrslice + group)
    grmed   <- tapply(sqres, group, quantile, probs=lts.quantile)
    meds    <- grmed[as.character(group)]
    stopifnot(!any(is.na(meds)))
    sel     <- (sqres <= meds)

    params[,,,lts.iter] <- pstart <- par

    ## Convergence check
    ## after a suggestion from David Kreil, kreil@ebi.ac.uk
    if(!is.null(cvg.check)) {
      cvgc    <- max(abs((hy - oldhy)/diff(range(hy))))
      cvgcCnt <- ifelse( cvgc < cvg.check$eps, cvgcCnt + 1, 0 )
      if (verbose)
        cat(sprintf("iter %2d: cvgc=%.5f%%, par=", as.integer(lts.iter), cvgc),
            sapply(par, function(x) sprintf("%9.3g",x)),"\n")
      if (cvgcCnt >= cvg.check$n)
        break
      oldhy <- hy
    }

  } ## end of for-loop (iter)
  if(verbose)
    cat("\n")

  ## Prepare the return result: an exprSet
  ## The transformed data goes into slot exprs.
  ## If input was allready an exprSet, keep the values of all the other slots.
  ## To the slot description@preprocessing, append the parameters and the
  ##    trimming selection.

  res <- descr <- NULL
  if (is(intensities, "exprSet")) {
    res <- intensities
    if (is(description(intensities), "MIAME")) {
      descr <- description(intensities)
    }
  }
  if(is.null(descr))   descr <- new("MIAME")
  if(is.null(res))     res   <- new("exprSet", description=descr)

  exprs(res) = hy[reord, ]
  if (describe.preprocessing) {
    vsnPreprocessing =  list(vsnParams        = par,
                           vsnParamsIter    = params,
                           vsnTrimSelection = sel[reord])
    class(vsnPreprocessing) = c("vsnPreprocessing", class(vsnPreprocessing))
    res@description@preprocessing = append(res@description@preprocessing, vsnPreprocessing)
  }
  return(res)
}
