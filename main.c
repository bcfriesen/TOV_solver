#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

struct param // parameter list to be passed to ODEs
{
  // pinit: central pressure
  // Gamma: polytropic index
  double Gamma, pinit;
  // single_star: if 0 (true), print grid along r
  int single_star;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  /* temperature. not part of ODE system being integrated so we have
     to leave it in the tag-along struct */
  double temp;
};

int make_grid (struct param params, double r, double r1, double y[]);

/* EOS function so we don't have to hard-code it into the structure of
   the ODEs. For now, since the 2 ODEs are for P and M, we'll assume
   the EOS is in the form rho = rho(P), which, I think, is usually the
   inverse of the structure of most real EOS routines. I'm doing it
   this way just so I don't have to change one of the ODE variables
   from P to rho. In the long run that would be the better way to do
   this, just so that it would play nicer with external EOSs. */
double eos_p_or_d(double p_or_d, struct param myparams, int d_to_p)
{
  return (pow(p_or_d, 1.0/myparams.Gamma)); // polytrope
}


// the system of ODEs to be integrated
int func (double r, const double y[], double f[], void *params)
{
  struct param myparams = *(struct param *)params;
  double (*p_or_d)(double, struct param, int); // EOS pointer
  p_or_d = &eos_p_or_d;

  // for a single star we can print P and M as functions of r
  if (myparams.single_star == 0)
  {
    printf("%12f %12f %12f\n", r, y[0], y[1]);
  }

  if (y[1] < 0.0) // can't have negative pressure
  {
    printf("%12f %12f %12f\n", r, y[0], myparams.pinit);
    return GSL_EBADFUNC; // this flag tells GSL integrator to quit
  } 

  // mass conservation equation
  f[0] = 4.0*M_PI*pow(r, 2.0)*p_or_d(y[1], myparams, 0);

  // TOV equation
  f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.0))
    * (p_or_d(y[1], myparams, 0)
    + (y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0)))
    + (y[0] + 4 * M_PI * pow(r, 3.0)
    * y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0))
    / (1.0 - (2.0 * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0]
    / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0) * r)));

  // ODEs tend to return NaNs at the same radius where P < 0
  if ((gsl_isnan(f[0]) || gsl_isnan(f[1]))) 
  {
    /* print values at the radius where ODEs diverge. this is almost
       always at the surface */
    printf("%12f %12f %12f\n", r, y[0], myparams.pinit);
    return GSL_EBADFUNC;
  }

  return GSL_SUCCESS;
}


int main (void)
{
  int i, status;
  struct param params;
  double r, r1, y[2];
  const double tiny = 1.0e-5, pmin = tiny, pmax = 0.2;
  const int MAX = 1000;
  // pressure/density arrays from tabulated EOS data
  double *eos_pres = NULL, *eos_dens = NULL;
  double *tmp1 = NULL, *tmp2 = NULL;
  int n_eos_pts;
  FILE *fp;

  fp = fopen("bck.eos", "r");

  if (fp == NULL) {
    fprintf(stderr, "Can't open file!\n");
    exit(1);
  }

  // read in # of EOS data points
  fscanf(fp, "%i", &n_eos_pts);

  printf("# of EOS data points: %i\n", n_eos_pts);
  eos_dens = (double *) malloc (n_eos_pts*sizeof(double));
  eos_pres = (double *) malloc (n_eos_pts*sizeof(double));
  tmp1 = (double *) malloc (n_eos_pts*sizeof(double));
  tmp2 = (double *) malloc (n_eos_pts*sizeof(double));
  for (i = 0; i < n_eos_pts; i++) {
    fscanf(fp, "%e %e %e %e",
	   &eos_dens[i], &eos_pres[i], &tmp1, &tmp2);
  }
  // read in density and pressure data

  fclose(fp);

  // set up interpolators for EOS
  //  params.acc = gsl_interp_accel_alloc ();
  //  params.spline = gsl_spline_alloc (gsl_interp_cspline, n_eos_pts);

  printf ("%12s %12s %12s\n", "R", "M(r=R)", "P(r=0)");

  r = 1.0e-5; // the integrator will go nuts if we start right at r=0
  r1 = 3.0; // some final 'radius' (much larger than actual radius)
  params.single_star = 1;
  for (i = 0; i < MAX; i++)
  {
    /* central mass always starts at 0 (or very close, for numerical
       reasons) */
    y[0] = 1.0e-6;
    y[1] = pmin + (double)i*((pmax - pmin)/(double)MAX);
    params.pinit = y[1];
    /* This function is useful if you want to plot, e.g., central
     * pressure vs. total mass. You can also hang onto the run of
     * pressure with radius, which can be interesting when compared to
     * the Newtonian case. */
    status = make_grid(params, r, r1, y);
  }

  printf("\nnow for a single star!\n");
  printf("%12s %12s %12s\n", "r", "M(r)", "P(r)");
  // print P and M vs r for Chandrasekhar-mass star
  y[0] = 1.0e-6;
  y[1] = 0.03; // pinit (roughly) for maximum mass
  params.single_star = 0;
  params.pinit = y[1];
  status = make_grid(params, r, r1, y);

  //  gsl_spline_free (params.spline);
  //  gsl_interp_accel_free (params.acc);
  free (eos_pres);
  free (eos_dens);
  free (tmp1);
  free (tmp2);

  eos_dens = NULL;
  eos_pres = NULL;
  tmp1 = NULL;
  tmp2 = NULL;  
  return 0;
}


int make_grid (struct param params, double r, double r1, double y[])
{
  int status;
  int *fake_jac;
 /* The integrator RK8PD (Runge-Kutta Prince-Dormand) doesn't need the
  * jacobian (some more sophisticated integrators do) so make this a
  * null pointer.  The linker will complain that the pointer type is
  * incompatable since it doesn't point to a function with a bunch of
  * arguments, but the pointer won't be used anyway so it doesn't
  * matter. */
  fake_jac = 0;
		   
  gsl_odeiv2_system sys = {func, fake_jac, 2, &params};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new
                         (&sys, gsl_odeiv2_step_rk4, 1.0e-6, 1.0e-6, 0.0);

  status = gsl_odeiv2_driver_apply(d, &r, r1, y);

  gsl_odeiv2_driver_free (d);
  return 0;
}
