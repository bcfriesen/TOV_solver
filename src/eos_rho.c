#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include "eos_rho.h"
#include "root_func.h"
#include "global.h"

// "reverse" EOS. input: P, T. output: rho
double
eos_rho (double pres, void *params)
{
  struct param *eos_rho_params = (struct param *) params;
  int status;
  double rho_low, rho_hi, root;
  // root-finding stuff for inverting EOS
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s = NULL;
  gsl_function F;
  // set up root-finding machinery for inverting EOS
  F.function = &root_func;
  F.params = params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  status = gsl_root_fsolver_set (s, &F, rho_root_min, rho_root_max);

  tmp_pres = pres;
  do
    {
      status = gsl_root_fsolver_iterate (s);
      root = gsl_root_fsolver_root (s);
      rho_low = gsl_root_fsolver_x_lower (s);
      rho_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (rho_low, rho_hi, 0, 0.001);
    }
  while (status == GSL_CONTINUE);

  gsl_root_fsolver_free (s);
  return (root);
}
