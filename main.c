#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <math.h>

int func (double t, const double y[], double f[], void *params)
{
  double Gamma = *(double *)params;
  printf("r = %12f; y[0] = %12f; y[1] = %12f\n", t, y[0], y[1]);
  f[0] = 4.0*M_PI*pow(t, 2.0)*pow(y[1], 1.0/Gamma);
  f[1] = -(pow(y[1], 1.0/Gamma)*y[0]/pow(t, 2.0))*(1.0 + pow(y[1], (Gamma - 1.0)/Gamma))*(1.0 + (4.0*M_PI*y[1]*pow(t, 3.0))/y[0])*
         pow(1.0 - 2.0*y[0]/t, -1.0);
 // printf("f[0] = %12f; f[1] = %12f\n", f[0], f[1]);
  return GSL_SUCCESS;
}

int main (void)
{
  double Gamma = 5.0/3.0;
  int *fake_jac = NULL;
  double t = 1.0e-5, t1 = 2.0;
  double y[2] = { 1.0e-6, 1.0 };

  gsl_odeiv2_system sys = {func, fake_jac, 2, &Gamma};

  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

  int status = gsl_odeiv2_driver_apply(d, &t, t1, y);

  gsl_odeiv2_driver_free (d);
  return 0;
}
