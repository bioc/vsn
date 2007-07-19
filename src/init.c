#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "vsn.h"

/* Registration information for DLL */
static const R_CallMethodDef callEntries[] = {
    { "vsn_c", (DL_FUNC) &vsn_c, 4}, 
    { "vsn2_point", (DL_FUNC) &vsn2_point, 5},
    { "vsn2_optim", (DL_FUNC) &vsn2_optim, 6},
    { "vsn2_trsf",  (DL_FUNC) &vsn2_trsf, 3},
    { "vsn2_scalingFactorTransformation", (DL_FUNC) &vsn2_scalingFactorTransformation, 1},
    { NULL, NULL, 0 }
};

/* Fow now, I here only describe the entry points for .Call, 
   not for .C, .Fortran, .External */
void R_init_vsn(DllInfo *info) {
  R_registerRoutines(info, NULL, callEntries, NULL, NULL);  
  R_useDynamicSymbols(info, FALSE);
}

void R_unload_vsn(DllInfo *info) {
  /* Here could be code to release resources. */
}

