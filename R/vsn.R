##-----------------------------------------------------------------
## Robust calibration and variance stabilization
## (C) Wolfgang Huber 2002-2003
## w.huber@dkfz.de
##-----------------------------------------------------------------
##----------------------------------------------------------------------
## vsn: the main function of this library
## This is one big chunk of a function - but lots of screen space are
## comments, and if you have suggestions on how to improve this function
## by splitting it up into smaller pieces, please come forward!
##----------------------------------------------------------------------
vsn <- function(intensities,
               lts.quantile = 0.5,
               verbose      = TRUE,
               niter        = 10,
               cvg.check    = NULL,
               pstart       = NULL,
               describe.preprocessing = TRUE)
{
  y <- getIntensityMatrix(intensities, verbose)

  ## Bureaucracy: make sure are the arguments are valid and plausible
  mess <- badarg <- NULL
  if (!is.numeric(lts.quantile) || (length(lts.quantile)!=1) ||
     (lts.quantile<0.5) || (lts.quantile>1)) {
    badarg <- "lts.quantile"
    mess   <- "Please specify a scalar between 0.5 and 1."
  }
  if (!is.numeric(niter) || (length(niter)!=1) || (niter<1)) {
    badarg <- "niter"
    mess   <- "Please specify a number >=1."
  }
  if (!is.null(pstart)) {
    if (!is.numeric(pstart) || length(pstart) != 2*ncol(y) || any(is.na(pstart))) {
      badarg <- "pstart"
      mess   <- paste("Please specify a numeric vector of length", 2*ncol(y))
    }
    if (any(pstart[(ncol(y)+1):(2*ncol(y))]<=0)) {
      badarg <- "pstart"
      mess   <- "Please specify non-negative values for the factors."
    }
  }
  if (!is.logical(verbose)) {
    badarg <- "verbose"
    mess   <- "Please specify a logical value."
  }
  ## Error handling
  if(!is.null(mess))
    stop(paste("The argument ", badarg, " has an invalid value:", get(badarg),
               "\n", mess, "\n", sep=""))

  ## Print welcome message
  if (verbose)
    cat("vsn: ", nrow(y), " x ", ncol(y), " matrix. Please wait for ",
        niter+1, " dots:\n.", sep="")

  ##----------------------------------------------------------------------
  ## guess a parameter scale, set boundaries for optimization,
  ## and, if they are not user-supplied, set the start parameters
  ##----------------------------------------------------------------------
  pscale  = plower = numeric(2*ncol(y))
  pscale[1:ncol(y)] = 1
  for (j in 1:ncol(y)) {
    pscale[ncol(y)+j] = 1/diff(quantile(y[,j], probs=c(0.25, 0.75)))
  }
  ## lower boundary for 'factors': a small positive value
  ## no boundary for 'offsets'
  plower = c(rep(-Inf,ncol(y)), pscale[(ncol(y)+1):(2*ncol(y))]/1e8)
  if (is.null(pstart))
    pstart = c(rep(0,ncol(y)), pscale[(ncol(y)+1):(2*ncol(y))])

  ## factr controls the convergence of the "L-BFGS-B" method. Convergence
  ## occurs when the reduction in the objective is within this factor of
  ## the machine tolerance. Default is 1e7, that is a tolerance of about
  ## 1e-8. Here we use 5e8 to save a little time.
  control     = list(trace=0, maxit=4000, parscale=pscale, factr=5e8)
  optim.niter = 10

  ## a place to save the trajectory of estimated parameters along the iterations:
  params  = matrix(NA, nrow=length(pscale), ncol=niter)
  rownames(params) = c(paste("offs", 1:ncol(y), sep=""), paste("fac", 1:ncol(y), sep=""))
    
  ## Workspace. This has two purposes: 1. we do not need to pass y around,
  ## 2. the intermediate results of "ll" can be used by "grll" (see functions
  ## vsnll, vsngrll)
  ws = new.env(hash=TRUE)

  ##----------------------------------------------------------------------------
  ## In the following, we define two functions: ll, grll. Doing this inside the 
  ## function vsn is one simply way to hide them from outside
  ##----------------------------------------------------------------------------
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
  ## Return value:
  ## Profile Log-Likelihood
  ##
  ## Side effects:
  ## 1. The function expects the matrix of input data y_ki in the environment ws.
  ## Since the matrix can be large (O(10^4) x O(10^2)), we want to avoid passing
  ## around this matrix by value.
  ## 2. The function stores intermediate results in that environment, too.
  ## They can then be re-used by grll.
  ##
  ## The variables stored in environment ws are
  ## y      the matrix of the y_ki
  ## nry    number of rows of y
  ## ncy    number of columns of y
  ## ly     the matrix of the a_i + b_i * y_ki
  ## asly   the matrix of the asinh(a_i + b_i * y_ki)
  ## res    the matrix of the residuals epsilon_ki
  ## ssq    sum of squared residuals sum(res^2)
  ## p      the call arguments of ll() . With this, when grll()
  ##        is called it can double-check whether it is indeed called
  ##        with the same arguments (see below for details)
  ##----------------------------------------------------------
  ll = function(p) {
    assign("p", p, envir=ws)
    with(ws, {
      offs = matrix(p[     1:ncy ], ncol=ncy, nrow=nry, byrow=TRUE)
      facs = matrix(p[ncy+(1:ncy)], ncol=ncy, nrow=nry, byrow=TRUE)
      ly   = offs + facs * y
      asly = asinh(ly)
      res  = asly - rowMeans(asly) ## residuals
      ssq  = sum(res*res)
      rv   = nry*ncy/2*log(ssq) - sum(log(facs/sqrt(1+ly*ly)))
    } )
    return(get("rv", ws))
  }
  ##--------------------------------------------------------------
  ## Gradient of the profile log likelihood
  ##--------------------------------------------------------------
  grll = function(p) {
    ## Generally, optim() will call the gradient of the objective function (gr)
    ## immediately after a call to the objective (fn) at the same parameter values.
    ## Anyway, we like to doublecheck
    if(any(p!=get("p", ws)))
      stop(paste(
        "\n\n\The function grll (likelihood-gradient) was called with different\n",
        "parameters than the previous call of ll (likelihood-value).\n",
        "This should never happen. Please contact the package maintainer:\n",
        "w.huber@dkfz.de\n\n"))

    with(ws, {
      dhda      = 1/sqrt(1+ly*ly)
      dlndhdyda = -ly/(1+ly*ly)
      gra       = nry*ncy/ssq*res*dhda - dlndhdyda
      rv        = c(colSums(gra), colSums(gra*y) - nry/p[(ncy+1):(ncy*2)])
    } )
    return(get("rv", ws))
  }
  
  sel     = rep(TRUE, nrow(y))   ## trimmed data points (LTS regression)
  oldh    = Inf  ## for calculating a convergence criterion: earlier result
  cvgcCnt = 0    ## counts the number of iterations that have already met the convergence criterion
  
  ##--------------------------------------------------
  ## begin of the outer LL iteration loop
  ##--------------------------------------------------
  for(lts.iter in 1:niter) {
    assign("y",   y[sel,],            envir=ws)
    assign("nry", length(which(sel)), envir=ws)
    assign("ncy", ncol(y),            envir=ws)

    p0 = pstart
    for (optim.iter in 1:optim.niter) {
      o  = optim(par=p0, fn=ll, gr=grll, method="L-BFGS-B",
                  control=control, lower=plower)
      if (o$convergence==0) next

      if(o$convergence==52) {
        ## ABNORMAL_TERMINATION_IN_LNSRCH
        ## This seems to indicate that a stepwidth to go along the gradient could not be found,
        ## probably because the start point p0 was already right at the optimum. Hence, try
        ## again from a slightly different start point
        ## cat("lts.iter=", lts.iter, "optim.iter=", optim.iter, "pstart was", p0, "now trying ")
        p0 = p0 + runif(length(pstart), min=0, max=0.01) * pscale
        ## cat(p0, "\n")
      } else if(o$convergence==1) {
        ## This seems to indicate that the max. number of iterations has been exceeded. Try again
        ## with more
        ## cat("lts.iter=", lts.iter, "optim.iter=", optim.iter, "maxit was", control$maxit, "now trying ")
        control$maxit = control$maxit*2
        ## cat(control$maxit, "\n")
      } else {
        stop(paste("Likelihood optimization: the function optim() returned the value convergence=",
                   o$convergence, "\nPlease make sure your data is good.",
                   "If so, contact the package maintainer.\n", sep=""))
      }
    }

    if (o$convergence!=0)
      stop(paste("Likelihood optimization did not converge even after", optim.niter, "calls to optim().",
                   "\nPlease make sure your data is good. If the problem persists,",
                   "\nplease contact the package maintainer.\n"))
    if (any(o$par[(ncol(y)+1):(2*ncol(y))]<0))
      stop(paste("Likelihood optimization produced negative parameter estimates in spite of constraints.",
                 "\nPlease contact the package maintainer.\n"))
    if(verbose)
      cat(".")

    # ----------------------------------------
    # selection of points in a LTS fashion
    # 1. calculate residuals; cf. ll()
    # ----------------------------------------
    offs   = matrix(o$par[    1      :   ncol(y) ], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
    facs   = matrix(o$par[(ncol(y)+1):(2*ncol(y))], nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
    asly   = asinh(offs + facs * y)
    res    = asly - rowMeans(asly)
    rsqres = rowSums(res*res)
    hmean  = rowSums(asly)

    # 2. select those data points within lts.quantile; do this separately for
    # each of nrslice slices along hmean
    nrslice = 5
    group   = ceiling(rank(hmean)/length(hmean)*nrslice)
    grmed   = tapply(rsqres, group, quantile, probs=lts.quantile)
    sel     = (rsqres <= grmed[group])

    params[, lts.iter] = pstart = o$par
    h = vsnh(y, o$par)

    ## Convergence check
    ## after a suggestion from D Kreil 2003, UCgenetics@Kreil.Org
    if(!is.null(cvg.check)) {
      cvgc    = max(abs((h - oldh)/range(h)))
      cvgcCnt = ifelse( cvgc < cvg.check$eps, cvgcCnt + 1, 0 )
      if (verbose)
        cat(sprintf("iter %2d: cvgc=%.5f%%, par=", as.integer(lts.iter), cvgc),
            sapply(o$par, function(x) sprintf("%9.3g",x)),"\n")
      if (cvgcCnt >= cvg.check$n)
        break
      oldh = h
    }
    
  } ## end of for-loop (iter)
  if(verbose)
    cat("\n")

  ## Prepare the return result: an exprSet
  ## The transformed data goes into slot exprs.
  ## If input was allready an exprSet, pass on the values all the other slots.
  ## To the slot description@preprocessing, append the parameters and the
  ##    trimming selection.  
  res = descr = NULL
  if (class(intensities)=="exprSet") {
    res = intensities
    if (class(description(intensities))=="MIAME") {
      descr = description(intensities)
    }
  }
  if(is.null(descr))   descr = new("MIAME")
  if(is.null(res))     res   = new("exprSet", description=descr)

  exprs(res) = vsnh(y, o$par)
  if (describe.preprocessing)
    res@description@preprocessing = append(res@description@preprocessing,
                      list(vsnParams        = params[,ncol(params)],
                           vsnParamsIter    = params,
                           vsnTrimSelection = sel))
  return(res)
}

