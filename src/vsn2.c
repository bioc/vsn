/*********************************************************************
  C functions for vsn 2.X
  (C) W. Huber 2002-2007
  This code replaces the deprecated code in vsn.c
  See the vignette 'incremental.Rnw' for the derivation of the maths
***********************************************************************/
#include <R.h>
#include <Rdefines.h>

#include <R_ext/Applic.h>         /* for lbfgsb */
#include <R_ext/Utils.h>          /* for R_CheckUserInterrupt */
extern double asinh(double);

#undef VSN_DEBUG 

typedef struct {
  double *y;       /* expression matrix: y_ik     */
  int nrow;        /* no. of features             */
  int ncol;        /* no. of chips                */
  int ntot;        /* no. of data points that are not NA (if none, this should be nrow*ncol) */
  int npar;        /* no. of parameters */

  int *strat;      /* strat[j] is the index of the first element of j-th stratum */
  int nrstrat;     /* no. of strata */

  int profiling;   /* 0 for normal likelihood, 0xffff for profiling                  */
  double *mu;      /* mu and sigma^2. For the normal likelihood, these are provided. */
  double sigsq;    /* For profile likelihood, they are computed from the data.       */

  /* Workspaces -  used to store intermediate results from the computation of the */
  /* likelihood function. These are reused in the computation of the gradient.    */
  double *ly;      /* is called Y_ki in the vignette   */
  double *asly;    /* is called h_ki in the vignette   */
  double *resid;   /* r_ki  */
  double *ma;      /* A_ki  */
  double *mb;      /* B_ki  */

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
    double z, a, b;
 
    nc = px->ncol;
    nr = px->nrow;
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
          hy[i + j*nr] = asinh((z+a)/b);
	}
      }
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
  double aj, bj, mu, s, z, ll, ssq, jac1, jac2, jacobian, scale, residuals;
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
    bj = b[j];
    ni = 0;
    for(i = px->strat[j]; i < px->strat[j+1]; i++){
      z = px->y[i];
      if(!ISNA(z)) {
	z = (z+aj)/bj;  
	px->ly[i]   = z;
	px->asly[i] = asinh(z);

	s = 1.0/sqrt(1.0+z*z); 
	px->ma[i] = s;
	px->mb[i] = z*s;
	
	jac1 += log(1.0+z*z); 
        ni++;
      } else {
	px->ly[i] = px->asly[i]	= px->ma[i] = px->mb[i] = NA_REAL;
      }
    } /* for i */
    jac2 += ni*log(bj);
    nt += ni;
  } /* for j */
  
  jacobian = jac1/2.0 + jac2;

  if(px->ntot != nt)   /* Just double-check - the code with 'nt' can be removed in future versions */
    error("Internal error 1 in 'loglik'.");

  /*---------------------------------------------------------------*/
  /* 2nd sweep through the data: compute r_ki                      */
  /*---------------------------------------------------------------*/
  ssq = 0.0;
  nt = 0;
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
        nt++;
      } else {
        z = NA_REAL;
      }
      px->resid[j*nr+i] = z;
    } /* for j */
  } /* for i */

  if(px->ntot != nt)   /* Just double-check - the code with 'nt' can be removed in future versions */
    error("Internal error 2 in 'loglik'.");

  if(px->profiling) {
    /* Negative profile log likelihood */
    residuals = (double)nt/2.0;
    scale = (double)nt/2.0 * log(2.0*M_PI*ssq/(double)nt);
    /* save for reuse in grad_loglik */
    px->sigsq = ssq/(double)nt;
  } else {
    /* Negative log likelihood */
    residuals = ssq /(2.0*(px->sigsq));
    scale = ((double)nt/2.0)*log(2.0*M_PI*(px->sigsq)); 
  }

  ll = scale + residuals + jacobian;

  /* Rprintf("loglik[%d]=%8g: scale=%8g, res=%8g, jac1=%8g, jac2=%8g, nt=%d, sigsq=%8g, ssq/nt=%8g\n", 
     px->profiling, ll, scale, residuals, jac1, jac2, nt, px->sigsq, ssq/(double)nt); */

  return(ll);
}

/*------------------------------------------------------------
   gradient
------------------------------------------------------------*/
void grad_loglik(int n, double *par, double *gr, void *ex)
{
  double *b;      
  vsn_data *px;
  double sa, sb, rfac, bki, z;
  int j, k, nj;

  px = (vsn_data*) ex;
  b  = par + px->nrstrat;

  for(j=0; j < px->npar; j++) {
    if (px->lastpar[j] != par[j]) {
      Rprintf("%d\t%g\t%g\n", j, px->lastpar[j], par[j]);
      error("Parameters in 'grad_loglik' are different from those in 'loglik'.");
    }
  }

  rfac = 1.0/(px->sigsq);
  
  /* Loop over the data */
  for(j = 0; j < px->nrstrat; j++) {
    sa = sb = 0.0;
    nj = 0;
    for(k = px->strat[j]; k < px->strat[j+1]; k++) {
      z = px->resid[k];
      if(!ISNA(z)){
        bki = px->mb[k];
	z = z*rfac + bki;
        sa += z*(px->ma[k]);
        sb += z*bki;
	nj++;
      }
    } /* for k */
    gr[j]               =  sa / b[j];
    gr[j+(px->nrstrat)] = ((double)nj - sb) / b[j];
  }


#ifdef VSN_DEBUG
  Rprintf("grad_loglik    "); 
  for(j=0; j < px->npar; j++) Rprintf(" %9g", gr[j]); 
  Rprintf("\n"); 
#endif 
  return;
}

/* This setup function is used by vsn2_optim, vsn2_point, vsn2_trsf 
   It checks the arguments Sy, Spar Strat and copies relevant bits 
   into (*px) */

void setupEverybody(SEXP Sy, SEXP Spar, SEXP Sstrat, vsn_data* px) 
{ 
  int i, nr, nc, nt;  
  double* y;
  SEXP dimy;

  if(fabs(asinh(1.5)-1.1947632172871)>1e-10)
    error("Your 'asinh' function does not seem to work right.");

  PROTECT(dimy = getAttrib(Sy, R_DimSymbol));
  if((!isReal(Sy)) | isNull(dimy) | (LENGTH(dimy)!=2))
    error("Invalid argument 'Sy', must be a real matrix."); 
  if(!isReal(Spar))
    error("Invalid argument 'Spar', must be a real vector.");
  if(!isInteger(Sstrat)) 
    error("Invalid argument 'Sstrat', must be integer.");

  px->npar    = LENGTH(Spar);
  px->strat   = INTEGER(Sstrat);

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
  px->mb      = (double *) R_alloc(nr*nc,  sizeof(double));
  px->lastpar = (double *) R_alloc(np, sizeof(double));

  /* parameters  */
  cpar = (double *) R_alloc(np, sizeof(double));
  for(i=0; i < 2*ns; i++) 
    cpar[i] = REAL(Spar)[i];
  for(i=ns; i < 2*ns; i++)
    if(cpar[i] <=0 )
      error("'Spar': factors must be >0.");
  
  return(cpar);
} 

/*------------------------------------------------------------
   vsn2_point: R interface for calculation of loglikelihood and gradient
   This is used by vsnLikelihood
------------------------------------------------------------*/
SEXP vsn2_point(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Smu, SEXP Ssigsq)
{
  double* cpar;
  SEXP res;
  vsn_data x;

  setupEverybody(Sy, Spar, Sstrat, &x);
  cpar = setupLikelihoodstuff(Sy, Spar, Sstrat, Smu, Ssigsq, &x);
 
  res = allocVector(REALSXP, x.npar+1);

  REAL(res)[0] = loglik(x.npar, cpar, (void*) &x);
  grad_loglik(x.npar, cpar, REAL(res)+1, (void*) &x);

  return(res);
}

/*------------------------------------------------------------
   vsn2_optim
------------------------------------------------------------*/
SEXP vsn2_optim(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Smu, 
                SEXP Ssigsq, SEXP Soptimpar)
{
  int i, lmm, fail, fncount, grcount, maxit, trace, nREPORT;
  int *nbd;
  double *cpar, *lower, *upper, *scale, *z;
  double factr, pgtol, fmin, low;
  char msg[60];

  SEXP res, namesres, vfail, coef, dimcoef, mu, sigsq;
  vsn_data x;

  lmm      = 20;   
  fail     = 0;

  if(!(isNewList(Soptimpar)&&(LENGTH(Soptimpar)==7)))
    error("Invalid argument: 'Soptimpar' must be a real vector of length 7.");

  factr    = REAL(getListElement(Soptimpar, "factr"))[0];
  pgtol    = REAL(getListElement(Soptimpar, "pgtol"))[0];
  low      = REAL(getListElement(Soptimpar, "lower"))[0];
  maxit    = INTEGER(getListElement(Soptimpar, "maxit"))[0]; 
  trace    = INTEGER(getListElement(Soptimpar, "trace"))[0];
 
  fncount  = 0;
  grcount  = 0;
  nREPORT  = 1;

  setupEverybody(Sy, Spar, Sstrat, &x);
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
    lower[i] = low;
    upper[i] = 0.;
    scale[i] = 1.;
  } 
  for(i=0; i<x.nrstrat; i++) {
    nbd[i]           = 0;
    nbd[i+x.nrstrat] = 1; /* lower bound for factors */
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
  INTEGER(dimcoef)[0] = x.npar/(x.ncol*2);
  INTEGER(dimcoef)[1] = x.ncol;
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
SEXP vsn2_trsf(SEXP Sy, SEXP Spar, SEXP Sstrat)
{
  int maxs, i; 
  SEXP res, dimres, dimy;
  vsn_data x;

  setupEverybody(Sy, Spar, Sstrat, &x);
 
  if (LENGTH(Sstrat) != x.nrow) 
    error("Length of 'Sstrat' must be the same as the number of rows of 'Sy'.");

  maxs = x.npar / (x.ncol*2);
  for(i=0; i<LENGTH(Sstrat); i++)
    if(x.strat[i]<1 || x.strat[i]>maxs) {
      Rprintf("x.strat[%d]=%d but should be >=1 and <=%d\n", i, x.strat[i], maxs);
      error("Invalid argument 'Sstrat'.");
    }
  
  PROTECT(res = allocVector(REALSXP, x.nrow*x.ncol));
  dimres = allocVector(INTSXP, 2);
  INTEGER(dimres)[0] = x.nrow;
  INTEGER(dimres)[1] = x.ncol;
  setAttrib(res, R_DimSymbol, dimres);

  calctrsf(&x, REAL(Spar), REAL(res));

  UNPROTECT(1);
  return(res);
}



