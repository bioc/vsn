
/* vsn2.c: */
SEXP vsn2_optim(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Srefh, SEXP Srefsigma);
SEXP vsn2_point(SEXP Sy, SEXP Spar, SEXP Sstrat, SEXP Srefh, SEXP Srefsigma);
SEXP vsn2_trsf(SEXP Sy, SEXP Spar, SEXP Sstrat);
SEXP vsn2_scalingFactorTransformation(SEXP Sb);

/* vsn.c */
SEXP vsn_c(SEXP e_y, SEXP e_par, SEXP e_strat, SEXP e_what);
