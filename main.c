#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <math.h>

int func (double t, const double y[], double f[], void *params)
{
  if (y[1] < 0.0) {return GSL_EBADFUNC;} // can't have negative pressure
    
  double Gamma = *(double *)params;
//  printf("r = %12f; M[0] = %12f; P[1] = %12f\n", t, y[0], y[1]);
  f[0] = 4.0*M_PI*pow(t, 2.0)*pow(y[1], 1.0/Gamma);
  f[1] = -(pow(y[1], 1.0/Gamma)*y[0]/pow(t, 2.0))*(1.0 + pow(y[1], (Gamma - 1.0)/Gamma))*(1.0 + (4.0*M_PI*y[1]*pow(t, 3.0))/y[0])*
         pow(1.0 - 2.0*y[0]/t, -1.0);
  if ((gsl_isnan(f[0]) || gsl_isnan(f[1]))) {return GSL_EBADFUNC;}

  return GSL_SUCCESS;
}


int main (void)
{
  double Gamma = 5.0/3.0;
  int *fake_jac = NULL;
  double t, t1;
  double y[2], pinit;
  const double tiny = 1.0e-5;
  const double pmin = tiny, pmax = 0.3;
  int i;
  const int MAX = 1000;
  int status;

  printf("%12s %12s %12s\n", "RADIUS", "M(r = R)", "P(r = 0)");
  for (i = 0; i < MAX; i++)
  {
    t = 1.0e-5;
    t1 = 3.0;
    y[0] = 1.0e-6; // mass always starts at 0
    y[1] = pmin + (double)i*((pmax - pmin)/(double)MAX);
    pinit = y[1];
//    printf("M[0] = %12f; P[1] = %12f\n", y[0], y[1]);
    status = make_grid(Gamma, &t, t1, y);
    printf("%12f %12f %12f\n", t, y[0], pinit);
  }
  
  return 0;
}


int make_grid (double Gamma, double *t, double t1, double y[])
{
  int status;
  int *fake_jac = NULL;
  gsl_odeiv2_system sys = {func, fake_jac, 2, &Gamma};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

  status = gsl_odeiv2_driver_apply(d, t, t1, y);

  gsl_odeiv2_driver_free (d);
  return 0;
}
