#include <gsl/gsl_spline.h>
#include "eos_pres.h"
#include "global.h"

// "forward" EOS. input: rho, T. output: P
double
eos_pres (double rho, void *params)
{
  struct param *eos_pres_params = (struct param *) params;
  return (gsl_spline_eval (spline, rho, acc));
  // return (pow(rho, 1.0/params.Gamma)); // polytrope
}
