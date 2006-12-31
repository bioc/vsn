/*******************************************************************
   C functions for vsn 2.X
   (C) W. Huber 2002-2007
   This code replaces the deprecated code in vsn.c
   (and has evolved from that)
*******************************************************************/
#include <R.h>
#include <Rdefines.h>

#include <R_ext/Applic.h>         /* for lbfgsb */
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>          /* for R_CheckUserInterrupt */
extern double asinh(double);

/*#define VSN_DEBUG */
#undef VSN_DEBUG

typedef struct {
  int *strat;      /* For what=0 and 1, strat[j] is the index of the first element 
                      of j-th stratum, and the length of this array is nrstrat.
                      For what=2, strat[j] is the stratum of the j-th probe, and
                      the length of this array is nrow.                   */
  int nrstrat;     /* no. of strata                                       */

  double *y;       /* expression matrix: y_ik     */
  int nrow;        /* no. of features             */
  int ncol;        /* no. of chips                */

  double *ly;      /* affine transformed matrix: offs_ik + facs_ik * y_ik */
  double *asly;    /* transformed expression matrix: asinh(ly)            */
  double *rcasly;  /* row centered version of asly                        */
  double *dh;      /* another auxilliary array                            */

  double *lastpar;
  int npar;         /* no. of parameters */
  double ssq;
} vsn_data;


/*--------------------------------------------------
  Apply the transformation to the matrix px->y
  Note: in contrast in to previous versions, the result
   is returned on the glog-scale to basis 2, not e.
  However, the computations further below (vsnloglik etc.) 
  are all done using the asinh and natural log scale 
  (since the tried and tested legacy code works well)
  --------------------------------------------------*/
void vsn2trsf(vsn_data *px, double* par, double *hy)
{
    int i, j, ns, s, nr, nc;
    double fac, off, h0;      
    double oolog2;
 
    oolog2 = 1.0/log(2.0);
    nc = px->ncol;
    nr = px->nrow;
    ns = px->npar / (nc*2);
    
    for(i=0; i <nr; i++) {
	s  = (px->strat[i]) - 1;
	h0 = log2(2*par[s + ns*nc]); 
	for(j=0; j<nc; j++) {
	    off          = par[s + j*ns];
	    fac          = par[s + j*ns + ns*nc];
            /*Rprintf("i=%d j=%d off=%g fac=%g\n", i, j, off, fac);*/
	    hy[i + j*nr] = oolog2*asinh(px->y[i+ j*nr] * fac + off) - h0;
	}
    }
    return;
}

/*----------------------------------------------------------------------
  The function to be optimized:
  offs[j] for j=0,..,nrstrat-1 are the offsets, lambda(facs[j]) the factors. 
  Here, lambda(x) is a monotonous, continuous, differentiable and strictly 
  positive function.
  By running the minimization on the parameters transformed by lambda, 
  the constraint that the factors need to be >0 is automatically satisfied. 
  For the gradient (see below), note that: 
    d/dx f(lambda(x)) = f'(lambda(x))*lambda'(x)
-----------------------------------------------------------------------*/
double lambda(double x)    { return(x*x);}
double dlambdadx(double x) { return(2.0*x);}
double invlambda(double y) { return(sqrt(y));} 

/* double lambda(double x)    { if(x<1) {return(exp(x-1));} else {return(x);}}
double dlambdadx(double x) { if(x<1) {return(exp(x-1));} else {return(1.0);}}
double invlambda(double y) { if(y<1) {return(log(y)+1.0);} else {return(y);}} */

/* double lambda(double x)    {return(x);}
double dlambdadx(double x) {return(1.0);}
double invlambda(double y) {return(y);} */

/* double lambda(double x)    {return(exp(x));}
double dlambdadx(double x) {return(exp(x));}
double invlambda(double y) {return(log(y));} */

double vsnloglik(int n, double *par, void *ex)
{
  double *facs, *offs;      
  double fj, oj, s, z, res, jac;  
  int i, j;
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
      z           = (px->y[i])*fj + oj; 
      px->ly[i]   = z;
      px->asly[i] = asinh(z);
      s           = 1.0/sqrt(1.0+z*z);   /* Jacobi term dh/dy   */
      px->dh[i]   = s;
      jac        += log(s);    /* logarithm since log-likelihood */

      /** FIXME **/
      if(!R_finite(z)){  
	Rprintf("j=%d, i=%d, z=%g, s=%g, y[i]=%g, fj=%g, oj=%g\n", j, i, z, s, px->y[i], fj, oj);
	for(j=0; j < px->npar; j++) Rprintf(" %6g", par[j]); 
	Rprintf("\n"); 
	error("Zapperlot z");
      }

      if(!R_finite(jac)){  
	Rprintf("j=%d, i=%d, z=%g, s=%g, y[i]=%g, fj=%g, oj=%g\n", j, i, z, s, px->y[i], fj, oj);
	for(j=0; j < px->npar; j++) Rprintf(" %6g", par[j]); 
	Rprintf("\n"); 
	error("Zapperlot jac");
      }
    } /* for i */
    jac += (px->strat[j+1] - px->strat[j])*log(fj);
  } /* for j */
  
  /* calculate ssq of residuals and rcasly  */
  /* rcasly  = row-centered version of asly */
  /* ssq     = sum_k sum_i rcasly_ik^2      */
  px->ssq = 0.0;
  for(i=0; i<nr; i++){
    s = 0.0;
    for(j=0; j < nc; j++){
      s += px->asly[j*nr+i];
    }
    s /= nc;
    for(j=0; j < nc; j++){
      z = px->asly[j*nr+i] - s;
      px->rcasly[j*nr+i] = z;
      px->ssq += z*z;
    }
  }

  /* the negative profile log likelihood */
  res = nr*nc*log(px->ssq)/2. - jac;

#ifdef VSN_DEBUG
  Rprintf("vsnloglik %9g", res); 
  for(j=0; j < px->npar; j++) Rprintf(" %9g", par[j]); 
  Rprintf("\n"); 
#endif

  return(res);
}

/*------------------------------------------------------------
   gradient of vsnloglik
   (see p.206 in notebook)
------------------------------------------------------------*/
void vsnloglikgrad(int n, double *par, double *gr, void *ex)
{
  double *facs;       
  double s1, s2, s3, s4, z1, z2, z3, vorfak; 
  int i, j, k;
  int nr, nc;
  vsn_data *px;

  px   = (vsn_data*) ex;
  facs = par + px->nrstrat;
  nr   = px->nrow;
  nc   = px->ncol;

  for(i=0; i < px->npar; i++) {
    if (px->lastpar[i] != par[i]) {
      Rprintf("%d\t%g\t%g\n", i, px->lastpar[i], par[i]);
      error("Parameters in 'vsnloglikgrad' are different from those in 'vsnloglik'.");
    }
  }

  vorfak = nr*nc/(px->ssq);

  for(j = 0; j < px->nrstrat; j++) {
    s1 = s2 = s3 = s4 = 0;
    for(k = px->strat[j]; k < px->strat[j+1]; k++) {
      z1 = (px->rcasly[k]) * (px->dh[k]); 
      s1 += z1;             /* deviation term for offset  */
      s2 += z1 * px->y[k];  /* deviation term for factor  */

      z2 = px->ly[k];
      z3 = z2/(1+z2*z2);    /* jacobi term                */
      s3 += z3;             /* jacobi term for offset     */
      s4 += z3 * px->y[k];  /* jacobi term for factor     */
    }
    s4 -= (px->strat[j+1] - px->strat[j])/lambda(facs[j]);
    gr[j]             =  vorfak * s1 + s3;
    gr[px->nrstrat+j] = (vorfak * s2 + s4) * dlambdadx(facs[j]); /* chain rule */
  }

#ifdef VSN_DEBUG
  Rprintf("vsnloglikgrad      "); 
  for(j=0; j < px->npar; j++) Rprintf(" %9g", gr[j]); 
  Rprintf("\n"); 
#endif 
  return;
}

/*------------------------------------------------------------
   The interface to R
   what=0: do the full parameters optimization
   what=1: just calculate the transformation and gradients on 
           given parameters
   what=2: just calculate the transformation on given parameters
------------------------------------------------------------*/
SEXP vsn2_c(SEXP e_y, SEXP e_par, SEXP e_strat, SEXP e_what)
{
  int i, nc, nr, what;

  int     lmm      = 10;   
  int     fail     = 0;
  double  factr    = 5e+7;  /* see below */
  double  pgtol    = 0; 
  int     fncount  = 0;
  int     grcount  = 0;
  int     maxit    = 40000;
  int     trace    = 0; /* 6; */
  int     nREPORT  = 1;
  double  fmin;
  char    msg[60];

  double  *cpar;
  double  *lower;
  double  *upper;
  double  *scale;
  int     *nbd;

  SEXP res, dimres, dimy;
  vsn_data x;

  /* check input arguments */
  PROTECT(dimy = getAttrib(e_y, R_DimSymbol));
 
  if((!isReal(e_y)) | isNull(dimy) | (LENGTH(dimy)!=2))
    error("Invalid argument 'e_y', must be a real matrix."); 
  if(!isReal(e_par))
    error("Invalid argument 'e_par', must be a real vector.");
  if(!isInteger(e_strat))
    error("Invalid argument 'e_strat', must be integer.");
  if(!isInteger(e_what) || LENGTH(e_what)!=1)
    error("Invalid argument 'e_what', must be integer of length 1.");
  what = INTEGER(e_what)[0];
  if(what<0 || what>2)
    error("Invalid argument 'e_what', must be 0, 1, or 2.");

  if(fabs(asinh(1.5)-1.1947632172871)>1e-10)
    error("Your 'asinh' function does not seem to work right.");

  /* assign length information and pointers to data areas into local workspace 'x' */
  x.npar    = LENGTH(e_par);
  x.strat   = INTEGER(e_strat);
  nr        = INTEGER(dimy)[0];
  nc        = INTEGER(dimy)[1];
  x.nrow    = nr;
  x.ncol    = nc;
  x.y       = REAL(e_y);

  /* workspaces for function and gradient calculation */
  x.ly      = (double *) R_alloc(nr*nc,  sizeof(double)); 
  if (what<=1) {
      x.nrstrat = LENGTH(e_strat)-1;
      if (2*x.nrstrat != x.npar) 
	  error("Unconformable size of arguments 'e_par', 'e_strat'.");
      if(x.strat[0] != 0)
	  error("First element of argument 'e_strat' must be 0.");
      if(x.strat[x.nrstrat] != nr*nc)
	  error("Last element of argument 'e_strat' must be equal to length of 'n_y'.");
      for(i=0; i<x.nrstrat; i++) {
	  if(x.strat[i+1]<= x.strat[i])
	      error("Elements of argument 'e_strat' must be in ascending order.");
      }
      /* if what is 0, the return value is a vector with 
         the parameters and with "ifail". If it is 1, it is 
         a vector with the value of the likelihood function 
         and its gradient. In addition, we need several
         workspace arrays for intermediate results */
      x.asly    = (double *) R_alloc(nr*nc,  sizeof(double));
      x.rcasly  = (double *) R_alloc(nr*nc,  sizeof(double));
      x.dh      = (double *) R_alloc(nr*nc,  sizeof(double));
      x.lastpar = (double *) R_alloc(x.npar, sizeof(double));
      res       = allocVector(REALSXP, x.npar+1);
  } else {
      x.nrstrat = x.npar/2;
      if (LENGTH(e_strat) != x.nrow) 
	  error("Length of 'e_strat' must be the same as the number of rows of 'e_y'.");
      lmm = x.npar / (nc*2);
      for(i=0; i<LENGTH(e_strat); i++)
	  if(x.strat[i]<1 || x.strat[i]>lmm) {
	      Rprintf("x.strat[%d]=%d but should be >=1 and <=%d\n", i, x.strat[i], lmm);
	      error("Invalid argument 'e_strat'.");
	  }
      /* if what is 2, the return value is a matrix of the same size 
         as y */
      res                = allocVector(REALSXP, nr*nc);
      PROTECT(dimres     = allocVector(INTSXP, 2));
      INTEGER(dimres)[0] = nr;
      INTEGER(dimres)[1] = nc;
      setAttrib(res, R_DimSymbol, dimres);
      UNPROTECT(1);
  }
  PROTECT(res);

  /* parameters:
     transform the factors using "lambda" (see above) */
  cpar = (double *) R_alloc(x.npar, sizeof(double));
  for(i=0; i<x.nrstrat; i++) 
    cpar[i] = REAL(e_par)[i];
  for(i=x.nrstrat; i < 2*x.nrstrat; i++) {
    if(REAL(e_par)[i] <=0 )
      error("'e_par': factors must be >0.");
    cpar[i] = (what<=1) ? invlambda(REAL(e_par)[i]) : REAL(e_par)[i];
  }

  /* parameter bounds and scale */
  switch(what) {
      case 0:
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
	  lbfgsb(x.npar, lmm, cpar, lower, upper, nbd, &fmin, vsnloglik, vsnloglikgrad, &fail,
		 (void *) &x, factr, pgtol, &fncount, &grcount, maxit, msg,
		 trace, nREPORT); 
	  /* write new values in result */
	  for(i=0; i < x.nrstrat; i++) 
	      REAL(res)[i] = cpar[i];
	  for(i=x.nrstrat; i < 2*x.nrstrat; i++) 
	      REAL(res)[i] = lambda(cpar[i]);
	  REAL(res)[x.npar] = (double) fail;
	  break;
      case 1:
	  /* just calculate function values and gradients; this is mostly for 
	     debugging, see script testderiv.R in inst/scripts directory  */
	  REAL(res)[0] = vsnloglik(x.npar, cpar, (void*) &x);
	  vsnloglikgrad(x.npar, cpar, REAL(res)+1, (void*) &x);
	  break;
      case 2:
	  vsn2trsf(&x, cpar, REAL(res));
          break;
  }

  UNPROTECT(2); /* dimy, res */
  return(res);
}




/* Registration information for DLL */
static R_CallMethodDef callMethods[] = {
    { "vsn2_c", ( DL_FUNC ) &vsn2_c, 4,
        /* { REALSXP, REALSXP, INTSXP, INTSXP } */ },
    { NULL, NULL, 0 }
};

extern void R_init_vsn2(DllInfo *info);
extern void R_unload_vsn2(DllInfo *info);

void R_init_vsn2(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_vsn2(DllInfo *info) {
  /* Release resources. */
}




/* Here's an excerpt from the Readme file for Lbfgsb.2.1

          L-BFGS-B (version 2.1)    March 1997
by
    Ciyou Zhu                    or           Jorge Nocedal
    ciyou@ece.nwu.edu                        nocedal@ece.nwu.edu


	L-BFGS-B is written in FORTRAN 77, in double precision.  The
user is required to calculate the function value f and its gradient g.
In order to allow the user complete control over these computations,
reverse communication is used.  The routine setulb.f must be called
repeatedly under the control of the variable task.  The calling
statement of L-BFGS-B is

      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)


	Following is a description of all the parameters used in this call.

c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave   */


