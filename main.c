#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <math.h>

struct param // parameter list to be passed to ODEs
{
  double Gamma, pinit;
};

// EOS function so we don't have to hard-code it into the structure of the
// ODEs. For now, since the 2 ODEs are for P and M, we'll assume the EOS is
// in the form rho = rho(P), which, I think, is usually the inverse of the
// structure of most real EOS routines. I'm doing it this way just so I
// don't have to change one of the ODE variables from P to rho. In the long
// run that would be the better way to do this, just so that it would play
// nicer with external EOSs.
double eos_rho(double pres, double Gamma);


double eos_rho(double pres, double Gamma)
{
  return (pow(pres, 1.0/Gamma)); // polytrope
}

int func (double t, const double y[], double f[], void *params)
{
  struct param myparams = *(struct param *)params;
  double (*rho)(double, double); // EOS pointer
  rho = &eos_rho;

  if (y[1] < 0.0) // can't have negative pressure
  {
    printf("%12f %12f %12f\n", t, y[0], myparams.pinit);
    return GSL_EBADFUNC; // this flag tells GSL integrator to quit
  } 

  // ODEs written explicit in terms of polytropic EOS. this should change
  // because it's not flexible.


  f[0] = 4.0*M_PI*pow(t, 2.0)*pow(y[1], 1.0/myparams.Gamma);
  f[1] = -(pow(y[1], 1.0/myparams.Gamma)*y[0]/pow(t, 2.0))*(1.0 +
           pow(y[1], (myparams.Gamma - 1.0)/myparams.Gamma))*(1.0 +
           (4.0*M_PI*y[1]*pow(t, 3.0))/y[0])*
           pow(1.0 - 2.0*y[0]/t, -1.0);

// ODEs tend to return NaNs at the same radius where P < 0
  if ((gsl_isnan(f[0]) || gsl_isnan(f[1]))) 
  {
    printf("%12f %12f %12f\n", t, y[0], myparams.pinit);
    return GSL_EBADFUNC;
  }

  return GSL_SUCCESS;
}


int main (void)
{
  struct param params;
  params.Gamma = 5.0/3.0;
  double t, t1;
  double y[2];
  const double tiny = 1.0e-5;
  const double pmin = tiny, pmax = 0.2;
  int i;
  const int MAX = 1000;
  int status;

  for (i = 0; i < MAX; i++)
  {
    t = 1.0e-5;
    t1 = 3.0;
    y[0] = 1.0e-6; // central mass always starts at 0
    y[1] = pmin + (double)i*((pmax - pmin)/(double)MAX);
    params.pinit = y[1];
    status = make_grid(params, &t, t1, y);
  }
  
  return 0;
}

int make_grid (struct param params, double *t, double t1, double y[])
{
  int status;
  int *fake_jac = NULL; // the integrator RK8PD doesn't need the jacobian
  gsl_odeiv2_system sys = {func, fake_jac, 2, &params};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new
                         (&sys, gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

  status = gsl_odeiv2_driver_apply(d, t, t1, y);

  gsl_odeiv2_driver_free (d);
  return 0;
}
