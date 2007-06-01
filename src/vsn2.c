/*******************************************************************
   C functions for vsn 2.X
   (C) W. Huber 2002-2007
   This code replaces the deprecated code in vsn.c
   (and has evolved from that)
*******************************************************************/
#include <R.h>
#include <Rdefines.h>

#include <R_ext/Applic.h>         /* for lbfgsb */
#include <R_ext/Utils.h>          /* for R_CheckUserInterrupt */
extern double asinh(double);

/*#define VSN_DEBUG */
#undef VSN_DEBUG

typedef struct {
  double *y;       /* expression matrix: y_ik     */
  int nrow;        /* no. of features             */
  int ncol;        /* no. of chips                */
  int ntot;        /* no. of data points that are not NA (if none, this should be nrow*ncol) */
  int npar;        /* no. of parameters */

  int *strat;      /* strat[j] is the index of the first element of j-th stratum */
  int nrstrat;     /* no. of strata */

  double ssq;

  double *refh;    /* reference values mu and sigma^2 */
  double refsigsq; 

  /* Workspaces -  used to store intermediate results from the computation 
     of the likelihood function. These are reused in the computation of the 
     function gradient */
  double *ly;      /* affine transformed matrix: offs_ik + facs_ik * y_ik */
  double *asly;    /* transformed expression matrix: asinh(ly)            */
  double *resid;   /* row centered version of asly                        */
  double *dh;      /* another auxilliary array                            */

  double *lastpar;
} vsn_data;


/*--------------------------------------------------
  Apply the transformation to the matrix px->y
  Note: in contrast in to previous versions, the result
   is returned on the glog-scale to basis 2, not e.
  However, the likelihood computations further below
  are all done using the asinh and natural log scale 
  --------------------------------------------------*/
void calctrsf(vsn_data *px, double* par, double *hy)
{
    int i, j, ns, s, nr, nc;
    double z, fac, off;      
 
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
	  off = par[s + j*ns];
	  fac = par[s + j*ns + ns*nc];
          hy[i + j*nr] = asinh(z*fac + off);
	}
      }
    }
    return;
}

/*----------------------------------------------------------------------
  The function to be optimized:
  offs[j] for j=0,..,nrstrat-1 are the offsets, facs[j] the factors. 
  
  For the normal likelihood, see the vignette 'incremental.Rnw'
  For the profile likelihood, see the SAGMB 2003 paper, and my grey 
   notebook p. 206 for the gradient.
-----------------------------------------------------------------------*/
double loglik(int n, double *par, void *ex)
{
  double *facs, *offs;      
  double fj, oj, mu, s, z, res, ssq, jac;  
  int i, j, ni, nt;
  int nr, nc;
  vsn_data *px;

  R_CheckUserInterrupt();

  px   = (vsn_data*) ex;
  offs = par;
  facs = par + px->nrstrat;
  nr   = px->nrow;
  nc   = px->ncol;

  for(i=0; i < px->npar; i++) 
    px->lastpar[i] = par[i];  

  jac = 0.0;
  for(j=0; j < px->nrstrat; j++){
    fj = facs[j];
    oj = offs[j];
    for(i = px->strat[j]; i < px->strat[j+1]; i++){
      z = px->y[i];
      if(!ISNA(z)) {
	z = z*fj + oj;  /* Affine transformation */
	s = 1.0/sqrt(1.0+z*z);   /* Jacobi term dh/dy */
	px->ly[i] = z;
	px->asly[i] = asinh(z);
	px->dh[i] = s;
	jac += log(s*fj);    /* logarithm since log-likelihood */
      } else {
	px->ly[i] = px->asly[i]	= px->dh[i] = NA_REAL;
      }
    } /* for i */
  } /* for j */
  
  /* calculate ssq of residuals and resid  */
  /* resid  = row-centered version of asly */
  /* ssq     = sum_k sum_i resid_ik^2      */
  ssq = 0.0;
  nt = 0;
  for(i=0; i<nr; i++){

    if(px->refh==NULL) {
      /* profiling: mu = arithmetic mean */
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
      nt += ni;
    } else {
      /* use the given parameter mu */
      mu = px->refh[i];
    }

    for(j=0; j < nc; j++){
      z = px->asly[j*nr+i];
      if(!(ISNA(mu)||ISNA(z))){
	z = z - mu;
	px->resid[j*nr+i] = z;
	ssq += z*z;
      } else {
	px->resid[j*nr+i] = NA_REAL;
      }
    }
  } /* for i */

  px->ssq = ssq;

  if(px->refh==NULL) {
    if(px->ntot != nt)      
      error("Internal error in 'loglik'.");
    /* Negative profile log likelihood */
    res = (px->ntot)*log(ssq)/2.0 - jac;
  } else {
    /* Negative log likelihood */
    /* Omitting the constant term nr*log(sqrt(2.0*M_PI)*sigma) 
       which is not relevant for the parameter optimization */
    res = ssq/(2.0*px->refsigsq) - jac;
  }

#ifdef VSN_DEBUG
  Rprintf("loglik %9g", res); 
  for(j=0; j < px->npar; j++) Rprintf(" %9g", par[j]); 
  Rprintf("\n"); 
#endif

  return(res);
}

/*------------------------------------------------------------
   gradient
------------------------------------------------------------*/
void grad_loglik(int n, double *par, double *gr, void *ex)
{
  double *facs;       
  double s1, s2, s3, s4, z1, z2, z3;
  double vorfak;
  int i, j, k;
  int nr, nc, nj;
  vsn_data *px;

  px   = (vsn_data*) ex;
  facs = par + px->nrstrat;
  nr   = px->nrow;
  nc   = px->ncol;

  for(i=0; i < px->npar; i++) {
    if (px->lastpar[i] != par[i]) {
      Rprintf("%d\t%g\t%g\n", i, px->lastpar[i], par[i]);
      error("Parameters in 'grad_loglik' are different from those in 'loglik'.");
    }
  }

  if(px->refh==NULL) {
    /* Negative profile log likelihood */
    vorfak = (px->ntot)/(px->ssq);
  } else {
    /* Negative log likelihood */
   vorfak = 1.0/(px->refsigsq);
  } 

  for(j = 0; j < px->nrstrat; j++) {
    s1 = s2 = s3 = s4 = 0.0;
    nj = 0;
    for(k = px->strat[j]; k < px->strat[j+1]; k++) {
      z1 = px->resid[k];
      if(!ISNA(z1)){
	z1 *= (px->dh[k]); 
	s1 += z1;             /* deviation term for offset  */
	s2 += z1 * px->y[k];  /* deviation term for factor  */
	  
	z2 = px->ly[k];
	z3 = z2/(1+z2*z2);    
	s3 += z3;             /* jacobi term for offset     */
	s4 += z3 * px->y[k];  /* jacobi term for factor     */
	nj++;
      }
    } /* for k */
    s4 -= (double)nj / facs[j];
    gr[j]             =  vorfak * s1 + s3;
    gr[px->nrstrat+j] = (vorfak * s2 + s4);
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

/* This setup function is used by vsn2_optim, vsn2_point (but not by vsn2_trsf) */
/* It sets up workspaces and processes the Sstrat parameters  */

double* setupLikelihoodstuff(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Srefh, SEXP Srefsigma, vsn_data* px) 
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

  /* Process Srefh and Srefsigma; If Srefh has length 0, then we do
     the 2002 model. If it has length nr, normalize against reference */
  if(!(isReal(Srefh)&&(isReal(Srefsigma))))
    error("Invalid arguments: 'Srefh' and 'Srefsigma' must be real vectors.");
  if((LENGTH(Srefh)==nr)&&(LENGTH(Srefsigma)==1)) {
    px->refh = REAL(Srefh);
    px->refsigsq = REAL(Srefsigma)[0];
    px->refsigsq *= px->refsigsq;  /* square */
  } else {
    if(LENGTH(Srefh)==0) {
      px->refh = NULL;
      px->refsigsq = 0.0;
    } else {
      error("Invalid length of arguments 'Srefh', 'Srefsigma'.");
    }
  }

  /* workspaces for function and gradient calculation */
  px->ly      = (double *) R_alloc(nr*nc,  sizeof(double)); 
  px->asly    = (double *) R_alloc(nr*nc,  sizeof(double));
  px->resid   = (double *) R_alloc(nr*nc,  sizeof(double));
  px->dh      = (double *) R_alloc(nr*nc,  sizeof(double));
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
SEXP vsn2_point(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Srefh, SEXP Srefsigma)
{
  double* cpar;
  SEXP res;
  vsn_data x;

  setupEverybody(Sy, Spar, Sstrat, &x);
  cpar = setupLikelihoodstuff(Sy, Spar, Sstrat, Srefh, Srefsigma, &x);
 
  res = allocVector(REALSXP, x.npar+1);

  REAL(res)[0] = loglik(x.npar, cpar, (void*) &x);
  grad_loglik(x.npar, cpar, REAL(res)+1, (void*) &x);

  return(res);
}

/*------------------------------------------------------------
   vsn2_optim
------------------------------------------------------------*/
SEXP vsn2_optim(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Srefh, 
                SEXP Srefsigma, SEXP Soptimpar)
{
  int i, lmm, fail, fncount, grcount, maxit, trace, nREPORT;
  int *nbd;
  double *cpar, *lower, *upper, *scale, *optimpar;
  double factr, pgtol, fmin;
  char msg[60];
  SEXP res;
  vsn_data x;

  lmm      = 20;   
  fail     = 0;

  if(!(isReal(Soptimpar)&&(LENGTH(Soptimpar)==4)))
    error("Invalid argument: 'Soptimpar' must be a real vector of length 4.");
  optimpar = REAL(Soptimpar);

  /* L-BFGS-B uses these two termination criteria:
     1. (f_k - f_k+1) / max(|f_k|, |f_k+1|, 1) <= factr * epsmch
       where epsmch is the machine precision.
     2. |gradient| < pgtol
     I use a combination of these two criteria, hoped to be a good compromise
     between not terminating too early on long flat valleys and never terminating
     because of round-off fluctuations.

     See L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
     Optimization, C. Zhu, R.H. Byrd, P. Lu and J. Nocedal (1996) */

  factr    = optimpar[0];  /* 4e4  */
  pgtol    = optimpar[1];  /* 2e-5 */
  maxit    = lrint(optimpar[2]);  /* 20000 */
  trace    = lrint(optimpar[3]);  /* 6     */
 
  fncount  = 0;
  grcount  = 0;
  nREPORT  = 1;

  setupEverybody(Sy, Spar, Sstrat, &x);
  cpar = setupLikelihoodstuff(Sy, Spar, Sstrat, Srefh, Srefsigma, &x);
 
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
    lower[i] = pgtol;
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

  /* write new values in result */
  res = allocVector(REALSXP, x.npar+1);
  for(i=0; i < x.npar; i++) 
    REAL(res)[i] = cpar[i];
  REAL(res)[x.npar] = (double) fail;
	  
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



