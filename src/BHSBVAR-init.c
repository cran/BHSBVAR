#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BHSBVAR_MAIN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_post_A_optim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_nonc_t(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_t(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_t_n(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_t_p(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BHSBVAR_MAIN",         (DL_FUNC) &_BHSBVAR_MAIN,         23},
    {"_BHSBVAR_post_A_optim", (DL_FUNC) &_BHSBVAR_post_A_optim, 13},
    {"_BHSBVAR_prior_nonc_t", (DL_FUNC) &_BHSBVAR_prior_nonc_t,  5},
    {"_BHSBVAR_prior_t",      (DL_FUNC) &_BHSBVAR_prior_t,       4},
    {"_BHSBVAR_prior_t_n",    (DL_FUNC) &_BHSBVAR_prior_t_n,     4},
    {"_BHSBVAR_prior_t_p",    (DL_FUNC) &_BHSBVAR_prior_t_p,     4},
    {NULL, NULL, 0}
};

void R_init_BHSBVAR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
