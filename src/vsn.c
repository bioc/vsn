/*******************************************************************
   optimization of the likelihood for vsn (with strata) in C. 

   The function vsnc expects:
   - data matrix
   - start parameters
   - vector strata (includes the number of genes until the start of a block)

   it returns:
   - new parameters after calling lbfgsb

*******************************************************************/
#include <R.h>
#include <Rdefines.h>

/* #define VSN_DEBUG */
#undef VSN_DEBUG

/* #ifndef HAVE_ASINH
#error "Don't have the 'asinh' function"
#endif
#ifndef HAVE_DECL_ASINH
#error "'asinh' function is not declared"
#endif */

/* Global variables */
int *strat;          /* strat[j] is index of first element of j-th stratum  */
double *y;           /* expression matrix: y_ik                             */
double *ly;          /* affine transformed matrix: offs_ik + facs_ik * y_ik */
double *asly;        /* transformed expression matrix: asinh(ly)            */
double *rcasly;      /* row centered version of asly                        */
double *dh;          /* another auxilliary array                            */
double ssq;
double *lastpar;

long int nrow;    /* no. of features   */
long int ncol;    /* no. of chips      */
long int nstrat;  /* no. of strata     */
long int npar;    /* no. of parameters */

/* For lbfgsb */
typedef double optimfn(int n, double *par, void *ex);
typedef void optimgr(int n, double *par, double *gr, void *ex);

void lbfgsb(int n, int lmm, double *x, double *lower, double *upper,
	    int *nbd, double *Fmin, optimfn fn, optimgr gr, int *fail,
	    void *ex, double factr, double pgtol, int *fncount,
	    int *grcount, int maxit, char *msg, int trace, int nREPORT);


/* faster than pow(x, 2) */
double sqr(double x) { return(x*x); }

/*------------------------------------------------------------
  Function to be optimized: likelihood

  offs[j] for j=0,..,nstrat-1 are the offsets,
  exp(facs[j]) the factors. By running the minimization 
  on the logarithm of the factors rather than on the factors
  themselves, the constraint factors>0 is automatically 
  satisfied. For the gradient, note that: 
    d/dx f(exp(x)) = f'(exp(x))*exp(x)
------------------------------------------------------------*/
double optfn(int n, double *par, void *ex)
{
  double *facs;      
  double *offs;      
  double s, res, jac;  
  long int i, j;

  offs = par;
  facs = par+nstrat;

  for(i=0; i<npar; i++) 
    lastpar[i] = par[i];  

  jac = 0.;
  for(j=0; j<nstrat; j++){
    for(i=strat[j]; i<strat[j+1]; i++){
      ly[i]   = y[i]*exp(facs[j]) + offs[j]; 
      asly[i] = asinh(ly[i]); 
      dh[i]   = 1.0/sqrt(1.0+sqr(ly[i]));
      jac    += facs[j] + log(dh[i]);    /* Jacobi term */
    }
  }

  /* calculate ssq and residuals            */
  /* rcasly  = row-centered version of asly */
  /* ssq     = sum_k sum_i rcasly_ik^2      */
  ssq = 0.;
  for(i=0; i<nrow; i++){
    s = 0.;
    for(j=0; j<ncol; j++){
      s += asly[j*nrow+i];
    }
    s /= ncol;
    for(j=0; j<ncol; j++){
      rcasly[j*nrow+i] = asly[j*nrow+i] - s;
      ssq += sqr(rcasly[j*nrow+i]);
    }
  }

  /* the negative profile log likelihood */
  res = nrow*ncol*log(ssq)/2. - jac;

  #ifdef VSN_DEBUG
  Rprintf("optfn %g", res); for(j=0; j<npar; j++) Rprintf(" %g", par[j]); Rprintf("\n"); 
  #endif

  return(res);
}

/*------------------------------------------------------------
   gradient of optfn
   (see p.206 in notebook)
------------------------------------------------------------*/
void optgr(int n, double *par, double *gr, void *ex)
{
  double *facs;       
  double s1, s2, s3, s4, z1, z2; 
  long int i, j, k;

  facs = par+nstrat;

  for(i=0; i<npar; i++) {
    if (lastpar[i] != par[i]) {
      Rprintf("%d\t%g\t%g\n", i, lastpar[i], par[i]);
      error("Parameters in 'optgr' are different from those in 'optfn'.");
    }
  }

  for(j=0; j<nstrat; j++) {
    s1 = s2 = s3 = s4 = 0;
    for(k=strat[j]; k<strat[j+1]; k++) {
      z1 = rcasly[k]*dh[k];       /* deviation term   */
      z2 = ly[k]/(1+sqr(ly[k]));  /* jacobi term      */
	s1 += z1;       /* deviation term for offset  */
	s2 += z1*y[k];  /* deviation term for factor  */
	s3 += z2;       /* jacobi term for offset     */
	s4 += z2*y[k];  /* jacobi term for factor     */
    }
    s4 -= (strat[j+1]-strat[j])/exp(facs[j]);
    gr[j]        = (nrow*ncol/ssq * s1 + s3);
    gr[nstrat+j] = (nrow*ncol/ssq * s2 + s4) * exp(facs[j]); /* chain rule */
  }

  #ifdef VSN_DEBUG
  Rprintf("optgr"); for(j=0; j<npar; j++) Rprintf(" %g", gr[j]); Rprintf("\n"); 
  #endif
  return;
}

/*------------------------------------------------------------
   The interface to R
------------------------------------------------------------*/
SEXP vsnc(SEXP e_y, SEXP e_par, SEXP e_strat, SEXP e_doopt)
{
  int i;

  int     lmm        = 10;   
  int     fail       = 0;
  double  factr      = 5e+7;  /* see below */
  double  pgtol      = 0; 
  int     fncount    = 0;
  int     grcount    = 0;
  int     maxit      = 40000;
  int     trace      = 0; /* 6; */
  int     nREPORT    = 1;
  double  fmin;
  char    msg[60];

  double  *cpar;
  double  *lower;
  double  *upper;
  double  *scale;
  int     *nbd;

  SEXP res, dimy;
 
  /* check input arguments */
  PROTECT(dimy = getAttrib(e_y, R_DimSymbol));
 
  if((!isReal(e_y)) | isNull(dimy) | (LENGTH(dimy)!=2))
    error("Invalid argument 'e_y', must be a real matrix."); 
  if(!isReal(e_par))
    error("Invalid argument 'e_par', must be a real vector.");
  if(!isInteger(e_strat))
    error("Invalid argument 'e_strat', must be integer.");
  if(!isLogical(e_doopt))
    error("Invalid argument 'e_doopt', must be logical.");

  if(abs(asinh(1.5)-1.1947632172871)>1e-10)
    error("Your 'asinh' function does not seem to work right.");

  /* assign length information and pointers to data areas into local caches */
  npar   = LENGTH(e_par);
  nstrat = LENGTH(e_strat)-1;
  nrow   = INTEGER(dimy)[0];
  ncol   = INTEGER(dimy)[1];
  y      = REAL(e_y);
  strat  = INTEGER(e_strat);

  /* check again */
  if (2*nstrat != npar) 
    error("Unconformable size of arguments 'e_par', 'e_strat'.");
  if(strat[0] != 0)
    error("First element of argument 'e_strat' must be 0.");
  if(strat[nstrat] != nrow*ncol)
    error("Last element of argument 'e_strat' must be equal to length of 'n_y'.");
  for(i=0; i<nstrat; i++) {
    if(strat[i+1]<=strat[i])
      error("Elements of argument 'e_strat' must be in ascending order.");
  }

  PROTECT(res = NEW_NUMERIC(npar+1));

  /* workspaces for function and gradient calculation */
  ly       = (double *) R_alloc(nrow*ncol, sizeof(double)); 
  asly     = (double *) R_alloc(nrow*ncol, sizeof(double));
  rcasly   = (double *) R_alloc(nrow*ncol, sizeof(double));
  dh       = (double *) R_alloc(nrow*ncol, sizeof(double));

  /* parameter bounds and scale */
  lower   = (double *) R_alloc(npar, sizeof(double));
  upper   = (double *) R_alloc(npar, sizeof(double));
  scale   = (double *) R_alloc(npar, sizeof(double));
  nbd     = (int *)    R_alloc(npar, sizeof(int));
  lastpar = (double *) R_alloc(npar, sizeof(double));
  cpar    = (double *) R_alloc(npar, sizeof(double));

  for(i=0; i<npar; i++) {
    lower[i]  = 0.;
    upper[i]  = 0.;
    scale[i]  = 1.;
    nbd[i]    = 0;   /* see below in the Readme file */
  } 
  for(i=0; i<nstrat; i++) 
    cpar[i] = REAL(e_par)[i];

  /* transform to log scale - see also comments for optfn */
  for(i=nstrat; i<2*nstrat; i++) {
    if(REAL(e_par)[i] <=0 )
      error("'e_par': factors must be >0.");
    cpar[i] = log(REAL(e_par)[i]);  
  }

  if(asLogical(e_doopt)) {
    /* optimize (see below for documentation of the function arguments) */
    lbfgsb(npar, lmm, cpar, lower, upper, nbd, &fmin, optfn, optgr, &fail,
	 (void *) NULL, factr, pgtol, &fncount, &grcount, maxit, msg,
	 trace, nREPORT); 
    /* write new values in result */
    for(i=0; i<nstrat; i++) 
      REAL(res)[i] = cpar[i];
    for(i=nstrat; i<2*nstrat; i++) 
      REAL(res)[i] = exp(cpar[i]);
    REAL(res)[npar] = (double) fail;

  } else {
    /* just calculate function values and gradients; this is mostly for 
       debugging, see script testderiv.R in inst/scripts directory  */
    REAL(res)[0] = optfn(npar, cpar, (void*) NULL);
    optgr(npar, cpar, REAL(res)+1, (void*) NULL);
  }

  UNPROTECT(2); /* dimy, res */
  return(res);
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

