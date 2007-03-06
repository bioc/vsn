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

  int *strat;      /* For what=0 and 1, strat[j] is the index of the first element 
                      of j-th stratum, and the length of this array is nrstrat.
                      For what=2, strat[j] is the stratum of the j-th probe, and
                      the length of this array is nrow.                   */
  int nrstrat;     /* no. of strata                                       */

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

double lambda(double x)    { return(x*x);}
double dlambdadx(double x) { return(2.0*x);}
double invlambda(double y) { return(sqrt(y));}

/* double lambda(double x)    { if(x<1) {return(exp(x-1));} else {return(x);}}
double dlambdadx(double x) { if(x<1) {return(exp(x-1));} else {return(1.0);}}
double invlambda(double y) { if(y<1) {return(log(y)+1.0);} else {return(y);}} */

/* double lambda(double x)    {return(fabs(x));}
double dlambdadx(double x) {return(x>0 ? 1.0 : -1.0);}
double invlambda(double y) {return(y);} */

/* double lambda(double x)    {return(exp(x));}
double dlambdadx(double x) {return(exp(x));}
double invlambda(double y) {return(log(y));} */

/*----------------------------------------------------------------------
  The function to be optimized:
  offs[j] for j=0,..,nrstrat-1 are the offsets, lambda(facs[j]) the factors. 
  Here, lambda(x) is a monotonous, continuous, differentiable and strictly 
  positive function.
  By running the minimization on the parameters transformed by lambda, 
  the constraint that the factors need to be >0 is automatically satisfied. 
  For the gradient, note that: 
    d/dx f(lambda(x)) = f'(lambda(x))*lambda'(x)
  
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
    fj = lambda(facs[j]);
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
    s4 -= (double)nj / lambda(facs[j]);
    gr[j]             =  vorfak * s1 + s3;
    gr[px->nrstrat+j] = (vorfak * s2 + s4) * dlambdadx(facs[j]); /* chain rule */
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
/* It sets up workspaces, processes the Sstrat parameters,
   and does the parameter transformation (lambda)  */

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

  /* parameters: transform the factors using "lambda" (see above) */
  cpar = (double *) R_alloc(np, sizeof(double));
  for(i=0; i < ns; i++) 
    cpar[i] = REAL(Spar)[i];
  for(i=ns; i < 2*ns; i++) {
    if(REAL(Spar)[i] <=0 )
      error("'Spar': factors must be >0.");
    cpar[i] = invlambda(REAL(Spar)[i]);
  }
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
SEXP vsn2_optim(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Srefh, SEXP Srefsigma)
{
  int i, lmm, fail, fncount, grcount, maxit, trace, nREPORT;
  int *nbd;
  double *cpar, *lower, *upper, *scale;
  double factr, pgtol, fmin;
  char msg[60];
  SEXP res;
  vsn_data x;

  lmm      = 20;   
  fail     = 0;

  /* L-BFGS-B uses these two termination criteria:
     1. (f_k - f_k+1) / max(|f_k|, |f_k+1|, 1) <= factr * epsmch
       where epsmch is the machine precision.
     2. |gradient| < pgtol
     I use a combination of these two criteria, hoped to be a good compromise
     between not terminating too early on long flat valleys and never terminating
     because of round-off fluctuations.

     See L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
     Optimization, C. Zhu, R.H. Byrd, P. Lu and J. Nocedal (1996) */
  factr    = 1e4; 
  pgtol    = 1e-5;
  maxit    = 200000;
 
  fncount  = 0;
  grcount  = 0;
  trace    = 0;    /* 6; */
  nREPORT  = 1;

  setupEverybody(Sy, Spar, Sstrat, &x);
  cpar = setupLikelihoodstuff(Sy, Spar, Sstrat, Srefh, Srefsigma, &x);
 
  lower   = (double *) R_alloc(x.npar, sizeof(double));
  upper   = (double *) R_alloc(x.npar, sizeof(double));
  scale   = (double *) R_alloc(x.npar, sizeof(double));
  nbd     = (int *)    R_alloc(x.npar, sizeof(int));
	  
  for(i=0; i<x.npar; i++) {
    lower[i]  = 0.;
    upper[i]  = 0.;
    scale[i]  = 1.;
    nbd[i]    = 0;   /* see below in the Readme file */
  } 
  
  /* optimize (see below for documentation of the function arguments) */
  lbfgsb(x.npar, lmm, cpar, lower, upper, nbd, &fmin, 
         loglik, grad_loglik, &fail,
	 (void *) &x, factr, pgtol, &fncount, &grcount, maxit, msg,
	 trace, nREPORT); 

  /* write new values in result */
  res = allocVector(REALSXP, x.npar+1);
  for(i=0; i < x.nrstrat; i++) 
    REAL(res)[i] = cpar[i];
  for(i=x.nrstrat; i < 2*x.nrstrat; i++) 
    REAL(res)[i] = lambda(cpar[i]);
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



