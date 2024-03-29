/*********************************************************************
  C functions for vsn 2.X
  (C) W. Huber 2002-2009
  This code replaces the deprecated code in vsn.c
  See the vignette 'likelihoodcomputations.Rnw' for the 
  derivation of the maths
***********************************************************************/

#include <R.h>
#include <Rdefines.h>

#include <R_ext/Applic.h>         /* for lbfgsb */
#include <R_ext/Utils.h>          /* for R_CheckUserInterrupt */
extern double asinh(double);

/* Transformation of the parameter 'b': see vignette */
#define FUN(b)  exp((b)) 
#define DFDB(b) exp((b)) 
/* #define FUN(b)  ((b<=0.0) ? exp(b) : b+1.0)  */
/* #define DFDB(b) ((b<=0.0) ? exp(b) : 1.0) */

#undef VSN_DEBUG

typedef struct {
  double *y;       /* expression matrix: y_ik     */
  int nrow;        /* no. of features             */
  int ncol;        /* no. of chips                */
  int ntot;        /* no. of data points that are not NA (if none, this should be nrow*ncol) */
  int npar;        /* no. of parameters */

  int *strat;      /* This array is used to store the mapping from rows to strata. Note that
                      at different times, it is used to store either the rows-to-strata mapping or 
                      the inverse, the strata-to-rows mapping. Using i for rows and j for strata:
                      In calctrsf, strat[i] is the stratum j (an integer identifier) of the i-th row.
                      In loglik and grad_loglik, strat[j] is the row index i of the first element 
                      of stratum j */
  int nrstrat;     /* Only used in loglik and grad_loglik: number of strata */

  int profiling;   /* 0 for normal likelihood, 0xffff for profiling                  */
  double *mu;      /* mu and sigma^2. For the normal likelihood, these are provided. */
  double sigsq;    /* For profile likelihood, they are computed from the data.       */
  int calib;       /* 0 for 'affine', all else for 'none'. */ 

  /* Workspaces -  used to store intermediate results from the computation of the */
  /* likelihood function. These are reused in the computation of the gradient.    */
  double *ly;      /* is called Y_ki in the vignette   */
  double *asly;    /* is called h_ki in the vignette   */
  double *resid;   /* r_ki  */
  double *ma;      /* A_ki  */

  double *lastpar;
} vsn_data;


/* From 'Writing R extensions' section 5.7.6 */
SEXP getListElement(SEXP list, char *str) {
  SEXP elmt=R_NilValue, names=getAttrib(list, R_NamesSymbol);
  int i;
  for (i=0; i<length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/*--------------------------------------------------
  Apply the transformation to the matrix px->y
  --------------------------------------------------*/
void calctrsf(vsn_data *px, double* par, double *hy)
{
    int i, j, ns, s, nr, nc;
    double z, a, b, fb;
 
    nc = px->ncol;
    nr = px->nrow;

    switch(px->calib) 
      {
      case 0:
	/* affine calibration */
	ns = px->npar / (nc*2);
    
	for(i=0; i <nr; i++) {
	  s  = (px->strat[i]) - 1; 
	  for(j=0; j<nc; j++) {
	    z = px->y[i+ j*nr];
	    if(ISNA(z)){
	      hy[i + j*nr] = NA_REAL;
	    } else {
	      a = par[s + j*ns];
	      b = par[s + j*ns + ns*nc];
	      hy[i + j*nr] = asinh(FUN(b)*z+a);
	    }
	  }
	}
	break;
      default:
	/* no calibration */
	a = par[0];
	fb = FUN(par[1]);
	for(i=0; i <nr*nc; i++) {
	  z = px->y[i];
	  hy[i] = ISNA(z) ? NA_REAL : asinh(fb*z+a);
	}
	break;
      }
    return;
}

/*----------------------------------------------------------------------
  The function to be optimized:
  a[j] for j=0,..,nrstrat-1 are the offsets, b[j] the factors. 
-----------------------------------------------------------------------*/
double loglik(int n, double *par, void *ex)
{
  double *a, *b;      
  double aj, bj, mu, z, ll, ssq, sigsq, jac1, jac2, jacobian, scale, residuals;
  int i, j, ni, nt;
  int nr, nc;
  vsn_data *px;

  R_CheckUserInterrupt();

  px = (vsn_data*) ex;
  a  = par;
  b  = par + px->nrstrat;
  nr = px->nrow;
  nc = px->ncol;

  for(i=0; i < px->npar; i++) 
    px->lastpar[i] = par[i];  

  /*---------------------------------------------------------------*/
  /* 1st sweep through the data: compute Y_ki, h(y_ki), A_ki, B_ki */
  /*---------------------------------------------------------------*/
  jac1 = jac2 = 0.0;
  nt = 0;
  for(j=0; j < px->nrstrat; j++){
    aj = a[j];
    bj = FUN(b[j]);

    /* Double-check - this expression can be removed in future versions */
    if(bj<=0)
      error("Nonpositive factor bj=%g (b[%d]=%g).\n", bj, j, b[j]);

    ni = 0;
    for(i = px->strat[j]; i < px->strat[j+1]; i++){
      z = px->y[i];
      if(!ISNA(z)) {
	z = bj*z+aj;  
	px->ly[i]   = z;
	px->asly[i] = asinh(z);
	px->ma[i]   = 1.0/sqrt(1.0+z*z); 
	jac1 += log(1.0+z*z); 
        ni++;
      } else {
	px->ly[i] = px->asly[i]	= px->ma[i] = NA_REAL;
      }
    } /* for i */
    jac2 += ni*log(bj);
    nt += ni;
  } /* for j */
  
  jacobian = jac1*0.5 - jac2;

  if(px->ntot != nt)   /* Double-check - the code with 'nt' can be removed in future versions */
    error("Internal error in 'loglik'.");

  /*---------------------------------------------------------------*/
  /* 2nd sweep through the data: compute r_ki                      */
  /*---------------------------------------------------------------*/
  ssq = 0.0;
  for(i=0; i<nr; i++){
    /* first, need to compute mu (if profiling) or just take it as 
       given (if not profiling) */
    if(px->profiling) {
      mu = 0.0;
      ni = 0;   /* count the number of data points in this row i (excl. NA) */
      for(j=0; j < nc; j++){
	z = px->asly[j*nr+i];
	if(!ISNA(z)){
	  mu += z;
	  ni++;
	}
      }
      mu = (ni>0) ? (mu/ni) : NA_REAL;
      px->mu[i] = mu;
    } else {
      /* no profiling: use the given parameter mu */
      mu = px->mu[i];
    }

    for(j=0; j < nc; j++){
      z = px->asly[j*nr+i];
      if(!(ISNA(mu)||ISNA(z))) {
        z = z-mu;
	ssq += z*z;
      } else {
        z = NA_REAL;
      }
      px->resid[j*nr+i] = z;
    } /* for j */
  } /* for i */

  if(px->profiling) {
    /* Negative profile log likelihood */
    /* Calculate sigsq and save for reuse in grad_loglik */
    sigsq = ssq/(double)nt; 
    px->sigsq = sigsq;
    residuals = (double)nt/2.0;
  } else {
    /* Negative log likelihood */
    sigsq = px->sigsq;
    residuals = ssq/(2.0*sigsq);
  }
  scale = (double)nt/2.0 * log(2.0*M_PI*sigsq); 
  ll = scale + residuals + jacobian;

#ifdef VSN_DEBUG
  Rprintf("negloglik[%d]=%8g: scale=%8g, res=%8g, jac1=%8g, jac2=%8g, nt=%d, sigsq=%8g, ssq/nt=%8g\n", 
     px->profiling, ll, scale, residuals, jac1, jac2, nt, px->sigsq, ssq/(double)nt); 
#endif

  return(ll);
}

/*------------------------------------------------------------
   gradient
------------------------------------------------------------*/
void grad_loglik(int n, double *par, double *gr, void *ex)
{
  double *b;      
  vsn_data *px;
  double sa, sb, rfac, aki, yki, lyki, z;
  int j, k, nj;

  px = (vsn_data*) ex;
  b  = par + px->nrstrat;

  for(j=0; j < px->npar; j++)
    if (px->lastpar[j] != par[j])
      error("Parameters in 'grad_loglik' are different from those in 'loglik': px->lastpar[%d]=%g but par[%d]=%g.\n",
             j, px->lastpar[j], j, par[j]);

  rfac = 1.0/(px->sigsq);
  
  /* Loop over the data */
  for(j = 0; j < px->nrstrat; j++) {
    sa = sb = 0.0;
    nj = 0;
    for(k = px->strat[j]; k < px->strat[j+1]; k++) {
      z = px->resid[k];
      if(!ISNA(z)){
        aki  = px->ma[k];       /* A_ki in the vignette */
        yki  = px->y[k];        /* y_ki in the vignette */
        lyki = px->ly[k];       /* Y_ki in the vignette */
	z = z*rfac + aki*lyki;
        sa += z*aki;
        sb += z*aki*yki;
	nj++;
      }
    } /* for k */
    gr[j]               = sa;
    gr[j+(px->nrstrat)] = DFDB(b[j])*(sb-(double)nj/FUN(b[j]));
  }


#ifdef VSN_DEBUG
  Rprintf(" gradient[%d]\npar:  ", px->profiling); 
  for(j=0; j < px->npar; j++) Rprintf(" %8g", par[j]); 
  Rprintf("\ngrad: "); 
  for(j=0; j < px->npar; j++) Rprintf(" %8g", gr[j]); 
  Rprintf("\n"); 
#endif

  return;
}

/* This setup function is used by vsn2_optim, vsn2_point, vsn2_trsf 
   It checks the arguments Sy, Spar Strat and copies relevant bits 
   into (*px) */

void setupEverybody(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Scalib, vsn_data* px) 
{ 
  int i, nr, nc, nt;  
  double* y;
  SEXP dimy;

  if(fabs(asinh(1.5)-1.1947632172871)>1e-10)
    error("Your 'asinh' function does not seem to work right.");

  PROTECT(dimy = getAttrib(Sy, R_DimSymbol));
  if((!isReal(Sy)) || isNull(dimy) || (LENGTH(dimy)!=2))
    error("Invalid argument 'Sy', must be a real matrix."); 
  if(!isReal(Spar))
    error("Invalid argument 'Spar', must be a real vector.");
  if(!isInteger(Sstrat)) 
    error("Invalid argument 'Sstrat', must be integer.");
  if(!isInteger(Scalib) || (LENGTH(Scalib)!=1)) 
    error("Invalid argument 'Scalib', must be integer of length 1.");

  px->npar    = LENGTH(Spar);
  px->strat   = INTEGER(Sstrat);
  px->calib   = INTEGER(Scalib)[0];

  px->y    = y  = REAL(Sy);
  px->nrow = nr = INTEGER(dimy)[0];
  px->ncol = nc = INTEGER(dimy)[1];

  nt = 0;
  for(i=0; i<nr*nc; i++)
    if(!ISNA(y[i]))
      nt++;
  px->ntot = nt;  

  UNPROTECT(1); 
  return;
}

/*------------------------------------------------------------------------------*/
/* This setup function is used by vsn2_optim, vsn2_point (but not by vsn2_trsf) */
/* It sets up workspaces and processes the Sstrat parameters                    */
/*------------------------------------------------------------------------------*/

double* setupLikelihoodstuff(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Smu, SEXP Ssigsq, vsn_data* px) 
{  
  int i, nr, nc, np, ns;
  double* cpar;

  nr = px->nrow;
  nc = px->ncol;
  np = px->npar;
  ns = px->nrstrat = LENGTH(Sstrat)-1;

  if ((2*ns) != np) 
    error("Wrong size of arguments 'Spar', 'Sstrat'.");
  if(px->strat[0] != 0)
    error("First element of argument 'Sstrat' must be 0.");
  if(px->strat[ns] != nr*nc)
    error("Last element of argument 'Sstrat' must be equal to length of 'n_y'.");
  for(i=0; i<ns; i++) {
    if((px->strat[i+1]) <= (px->strat[i]))
      error("Elements of argument 'Sstrat' must be in ascending order.");
  }

  /* Process Smu and Ssigsq; If Smu has length 0, then we do
     the 2002 model. If it has length nr, normalize against reference */
  if(!(isReal(Smu)&&(isReal(Ssigsq))))
    error("Invalid arguments: 'Smu' and 'Ssigsq' must be real vectors.");
  if((LENGTH(Smu)==nr)&&(LENGTH(Ssigsq)==1)) {
    px->mu    = REAL(Smu);
    px->sigsq = REAL(Ssigsq)[0];
    px->profiling = 0;
  } else {
    if(LENGTH(Smu)==0) {
      px->mu = (double *) R_alloc(nr, sizeof(double));
      px->sigsq = NA_REAL;
      px->profiling = 0xffff;
    } else {
      error("Invalid length of arguments 'Smu', 'Ssigsq'.");
    }
  }

  /* workspaces for function and gradient calculation */
  px->ly      = (double *) R_alloc(nr*nc,  sizeof(double)); 
  px->asly    = (double *) R_alloc(nr*nc,  sizeof(double));
  px->resid   = (double *) R_alloc(nr*nc,  sizeof(double));
  px->ma      = (double *) R_alloc(nr*nc,  sizeof(double));
  px->lastpar = (double *) R_alloc(np, sizeof(double));

  /* parameters  */
  cpar = (double *) R_alloc(np, sizeof(double));
  for(i=0; i < 2*ns; i++) 
    cpar[i] = REAL(Spar)[i];
  
  return(cpar);
} 

/*------------------------------------------------------------
   vsn2_point: R interface for calculation of loglikelihood and gradient
   This is used by vsnLikelihood
------------------------------------------------------------*/
SEXP vsn2_point(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Smu, SEXP Ssigsq, SEXP Scalib)
{
  double* cpar;
  SEXP res;
  vsn_data x;

  setupEverybody(Sy, Spar, Sstrat, Scalib, &x);
  cpar = setupLikelihoodstuff(Sy, Spar, Sstrat, Smu, Ssigsq, &x);
 
  res = PROTECT(allocVector(REALSXP, x.npar+1));

  REAL(res)[0] = loglik(x.npar, cpar, (void*) &x);
  grad_loglik(x.npar, cpar, REAL(res)+1, (void*) &x);

  UNPROTECT(1);
  return(res);
}

/*------------------------------------------------------------
   vsn2_optim
------------------------------------------------------------*/
SEXP vsn2_optim(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Smu, 
                SEXP Ssigsq, SEXP Soptimpar, SEXP Scalib)
{
  int i, lmm, fail, fncount, grcount, maxit, trace, nREPORT;
  int *nbd;
  double *cpar, *lower, *upper, *scale, *z;
  double factr, pgtol, fmin;
  char msg[60];

  SEXP res, namesres, vfail, coef, dimcoef, mu, sigsq;
  vsn_data x;

  lmm      = 5;   
  fail     = 0;
  nREPORT  = 10;

  if(!(isNewList(Soptimpar)&&(LENGTH(Soptimpar)==6)))
    error("Invalid argument: 'Soptimpar' must be a list of length 6.");

  factr    = REAL(getListElement(Soptimpar, "factr"))[0];
  pgtol    = REAL(getListElement(Soptimpar, "pgtol"))[0];
  maxit    = INTEGER(getListElement(Soptimpar, "maxit"))[0]; 
  trace    = INTEGER(getListElement(Soptimpar, "trace"))[0];
 
  fncount  = 0;
  grcount  = 0;

  setupEverybody(Sy, Spar, Sstrat, Scalib, &x);
  cpar = setupLikelihoodstuff(Sy, Spar, Sstrat, Smu, Ssigsq, &x);
 
  lower   = (double *) R_alloc(x.npar, sizeof(double));
  upper   = (double *) R_alloc(x.npar, sizeof(double));
  scale   = (double *) R_alloc(x.npar, sizeof(double));
  nbd     = (int *)    R_alloc(x.npar, sizeof(int));
	  
  /* nbd is an integer array of dimension n.
     On entry nbd represents the type of bounds imposed on the
     variables, and must be specified as follows:
        nbd(i)=0 if x(i) is unbounded,
               1 if x(i) has only a lower bound,
               2 if x(i) has both lower and upper bounds, and
               3 if x(i) has only an upper bound. */

  for(i=0; i<x.npar; i++) {
    scale[i] = 1.0;
  } 
  for(i=0; i<x.nrstrat; i++) {
    lower[i] = R_NegInf;
    upper[i] = R_PosInf;
    nbd[i]   = 0;
    lower[i+x.nrstrat] = -100.0;   /* exp(-100) is a very small positive number */
    upper[i+x.nrstrat] = +100.0;   /* exp(+100) is a very large positive number */
    nbd[i+x.nrstrat]   = 2;
  }

  /* optimize (see below for documentation of the function arguments) */
  lbfgsb(x.npar, lmm, cpar, lower, upper, nbd, &fmin, 
         loglik, grad_loglik, &fail,
	 (void *) &x, factr, pgtol, &fncount, &grcount, maxit, msg,
	 trace, nREPORT); 

  /* optimisation result "fail" */
  PROTECT(vfail = allocVector(INTSXP, 1));
  INTEGER(vfail)[0] = fail;

  /* optimisation result "sigsq" */
  PROTECT(sigsq = allocVector(REALSXP, 1));
  REAL(sigsq)[0] = x.sigsq;

  /* optimisation result "mu" */
  PROTECT(mu = allocVector(REALSXP, x.nrow));
  z = REAL(mu);
  for(i=0; i < x.nrow; i++) 
    z[i] = x.mu[i];
	  
  /* optimisation result "coef" */
  PROTECT(coef = allocVector(REALSXP, x.npar));
  z = REAL(coef);
  for(i=0; i < x.npar; i++) 
    z[i] = cpar[i];

  PROTECT(dimcoef = allocVector(INTSXP, 3)); 
  switch(x.calib)
    {
    case 0:
      INTEGER(dimcoef)[0] = x.npar/(x.ncol*2);
      INTEGER(dimcoef)[1] = x.ncol;
      break;
    default:
      INTEGER(dimcoef)[0] = 1;
      INTEGER(dimcoef)[1] = 1;
      break;
    }
  INTEGER(dimcoef)[2] = 2;
  setAttrib(coef, R_DimSymbol, dimcoef);

  /* return value: a list with four elements: fail, coefficients, mu, sigsq */
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, vfail);
  SET_VECTOR_ELT(res, 1, coef);  /* coefficients */
  SET_VECTOR_ELT(res, 2, sigsq); 
  SET_VECTOR_ELT(res, 3, mu);    

  PROTECT(namesres = allocVector(STRSXP, 4));
  SET_STRING_ELT(namesres, 0, mkChar("fail"));
  SET_STRING_ELT(namesres, 1, mkChar("coefficients"));
  SET_STRING_ELT(namesres, 2, mkChar("sigsq"));
  SET_STRING_ELT(namesres, 3, mkChar("mu"));
  setAttrib(res, R_NamesSymbol, namesres);

  UNPROTECT(7); /* vfail, sigsq, mu, coef, dimcoef, res, namesres */
  return(res);
}

/*------------------------------------------------------------
   vsn2_trsf
------------------------------------------------------------*/
SEXP vsn2_trsf(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Scalib)
{
  SEXP res, dimres;
  vsn_data x;

  setupEverybody(Sy, Spar, Sstrat, Scalib, &x);
 
  if (LENGTH(Sstrat) != x.nrow) 
    error("Length of 'Sstrat' must be the same as the number of rows of 'Sy'.");

  PROTECT(res = allocVector(REALSXP, x.nrow*x.ncol));
  dimres = allocVector(INTSXP, 2);
  INTEGER(dimres)[0] = x.nrow;
  INTEGER(dimres)[1] = x.ncol;
  setAttrib(res, R_DimSymbol, dimres);

  calctrsf(&x, REAL(Spar), REAL(res));

  UNPROTECT(1);
  return(res);
}

/*------------------------------------------------------------
   vsn2_scalingFactorTransformation
------------------------------------------------------------*/
SEXP vsn2_scalingFactorTransformation(SEXP Sb)
{
  int i, n;
  double *b, *r; 
  SEXP res;
  
  if(!isReal(Sb))
    error("Invalid argument 'Sb', must be a real vector.");

  n = LENGTH(Sb);
  b = REAL(Sb);
  
  /* No PROTECT since we do not call back into R or otherwise allocate memory */
  res = allocVector(REALSXP, n); 
  r = REAL(res);

  for(i=0; i<n; i++)
    r[i] = FUN(b[i]);

  return(res);
}



