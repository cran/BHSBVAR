#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _BHSBVAR_fevd_estimates(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_hd_estimates(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_irf_estimates(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_log_likelihood_function(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_MAIN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_beta(SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_ibeta(SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_nonc_t(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_t(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_t_n(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_prior_t_p(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BHSBVAR_sum_log_prior_densities(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BHSBVAR_fevd_estimates", (DL_FUNC) &_BHSBVAR_fevd_estimates,  7},
    {"_BHSBVAR_hd_estimates",   (DL_FUNC) &_BHSBVAR_hd_estimates,    6},
    {"_BHSBVAR_irf_estimates",  (DL_FUNC) &_BHSBVAR_irf_estimates,   6},
    {"_BHSBVAR_log_likelihood_function", (DL_FUNC) &_BHSBVAR_log_likelihood_function,  6},
    {"_BHSBVAR_MAIN",                    (DL_FUNC) &_BHSBVAR_MAIN,                    18},
    {"_BHSBVAR_prior_beta",              (DL_FUNC) &_BHSBVAR_prior_beta,               3},
    {"_BHSBVAR_prior_ibeta",             (DL_FUNC) &_BHSBVAR_prior_ibeta,              3},
    {"_BHSBVAR_prior_nonc_t",            (DL_FUNC) &_BHSBVAR_prior_nonc_t,             5},
    {"_BHSBVAR_prior_t",                 (DL_FUNC) &_BHSBVAR_prior_t,                  4},
    {"_BHSBVAR_prior_t_n",               (DL_FUNC) &_BHSBVAR_prior_t_n,                4},
    {"_BHSBVAR_prior_t_p",               (DL_FUNC) &_BHSBVAR_prior_t_p,                4},
    {"_BHSBVAR_sum_log_prior_densities", (DL_FUNC) &_BHSBVAR_sum_log_prior_densities,  4},
    {NULL, NULL, 0}
};

void R_init_BHSBVAR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
