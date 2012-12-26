#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>

// bounds for root-finder to invert EOS
// const double rho_root_min = 1.01e+6, rho_root_max = 4.9e+15; // use for bck.eos
// const double rho_root_min = 7.81, rho_root_max = 3.0e+16; // use for eosC
// const double rho_root_min = 1.1e+6, rho_root_max = 9.5e+15; // use for timmes.eos
const double rho_root_min = 1.1e+6, rho_root_max = 9.9e+14;	// use for helmholtz.eos


// interpolation machinery
gsl_interp_accel *acc = NULL;
gsl_spline *spline = NULL;

/* temporary storage of pressure when calling root-finder to invert
   the EOS */
double tmp_pres;

struct param			// parameter list to be passed to ODEs
{
  // pinit: central pressure
  // Gamma: polytropic index
  // rho_init: central density
  double Gamma, pinit, rho_init;
  // single_star: if 0 (true), print grid along r
  int single_star;
  /* temperature. not part of ODE system being integrated so we have
     to leave it in the tag-along struct */
  double temp;
};


int make_grid (struct param params, double r, double r1, double y[]);
double root_func (double rho, void *params);

// "forward" EOS. input: rho, T. output: P
double
eos_pres (double rho, void *params)
{
  struct param *eos_pres_params = (struct param *) params;
  return (gsl_spline_eval (spline, rho, acc));
  // return (pow(rho, 1.0/params.Gamma)); // polytrope
}


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

double
root_func (double rho, void *params)
{
  struct param *root_func_params = (struct param *) params;
  double (*pres) (double, void *);
  pres = &eos_pres;
  return (tmp_pres - pres (rho, root_func_params));
}

// the system of ODEs to be integrated
int
func (double r, const double y[], double f[], void *params)
{
  struct param *myparams = (struct param *) params;
  double (*pres) (double, void *);	// EOS pressure pointer
  double (*rho) (double, void *);	// EOS density pointer
  pres = &eos_pres;
  rho = &eos_rho;

// set a minimum pressure cutoff. if we don't, the ODE solver will wobble all
// over the surface and crash, or if you make the error tolerance really strict
// it'll integrate forever
  if (y[1] < 1.0e-9 * myparams->pinit)
    {
      if (myparams->single_star == 0)
	{
	  printf ("surface reached!\n");
	  printf ("%12s %12s %12s\n", "R (km)", "M (M_sun)",
		  "rho(0) (g/cm^3)");
	}
      printf ("%12le %12le %12le\n", r / 1.0e+5,
	      y[0] / GSL_CONST_CGSM_SOLAR_MASS, myparams->rho_init);
      return GSL_EBADFUNC;	// this flag tells GSL integrator to quit
    }
  // for a single star we can print M, P and rho as functions of r

  if (myparams->single_star == 0)
    {
      printf ("%12le %12le %12le %12le\n", r / 1.0e+5,
	      y[0] / GSL_CONST_CGSM_SOLAR_MASS, y[1], rho (y[1], myparams));
    }

  // mass conservation equation
  f[0] = 4.0 * M_PI * pow (r, 2.0) * rho (y[1], myparams);

  // TOV equation
  f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow (r, 2.0))
    * (rho (y[1], myparams)
       + (y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0)))
    * (y[0] + 4 * M_PI * pow (r, 3.0)
       * y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0))
    / (1.0 - (2.0 * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0]
	      / (pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0) * r)));

  // ODEs tend to return NaNs at the same radius where P < 0
  if ((gsl_isnan (f[0]) || gsl_isnan (f[1])))
    {
      /* print values at the radius where ODEs diverge. this is almost
         always at the surface */
      printf ("%12f %12f %12f\n", r / 1.0e+5,
	      y[0] / GSL_CONST_CGSM_SOLAR_MASS, myparams->rho_init);
      return GSL_EBADFUNC;
    }

  return GSL_SUCCESS;
}


int
main (void)
{
  int i, status;
  struct param params;
  double r, r1, y[2];
  // bounds for making grid of m(R) vs. rho(0)
  const double rho_min = 1.0e+13, rho_max = 9.0e+14;
  const int MAX = 100;
  // pressure/density arrays from tabulated EOS data
  int n_eos_pts;
  FILE *fp = NULL;
  double (*pres) (double, void *);	// EOS pressure pointer
  double (*rho) (double, void *);	// EOS density pointer
  double tmp;
  // pointers to tabulated EOS data
  double *eos_tab_pres = NULL, *eos_tab_dens = NULL;

  /* this is the "third arm" of the EOS inversion routine. GSL calls
     the function whose root is being found several times during
     initialization and if we leave the pressure here uninitialized,
     it might be 0 or something crazy, which will be out of bounds of
     the "forward" EOS, and since the forward EOS uses interpolation,
     an out-of-bounds value of pressure will make it crash. */
  tmp_pres = 1.0e+25;

  pres = &eos_pres;
  rho = &eos_rho;

  //  fp = fopen("bck.eos", "r");
  //  fp = fopen("eosC", "r");
  //  fp = fopen("timmes.eos", "r");
  fp = fopen ("helmholtz.eos", "r");

  if (fp == NULL)
    {
      fprintf (stderr, "Can't open file!\n");
      exit (1);
    }

  // read in # of EOS data points
  fscanf (fp, "%i", &n_eos_pts);

  eos_tab_dens = (double *) malloc (n_eos_pts * sizeof (double));
  eos_tab_pres = (double *) malloc (n_eos_pts * sizeof (double));

  // read in density and pressure data
  for (i = 0; i < n_eos_pts; i++)
    {
      fscanf (fp, "%le %le", &eos_tab_dens[i], &eos_tab_pres[i]);
    }
  fclose (fp);

  // set up interpolation machinery for EOS
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, n_eos_pts);
  gsl_spline_init (spline, eos_tab_dens, eos_tab_pres, n_eos_pts);

  // now make grid
  r = 1.0;			// the integrator will go nuts if we start right at r=0
  r1 = 1.0e+10;			// some final 'radius' (much larger than actual radius)

  printf ("%12s %12s %12s\n", "R", "M(r=R)", "rho(r=0)");

  params.single_star = 1;
  for (i = 0; i < MAX; i++)
    {
      // rho(0) evenly spaced in log
      params.rho_init = log10 (rho_min) +
	(double) i *((log10 (rho_max) - log10 (rho_min)) / (double) MAX);
      params.rho_init = pow (10.0, params.rho_init);

      y[1] = pres (params.rho_init, &params);
      y[0] = (4.0 / 3.0) * M_PI * pow (r, 3.0) * rho (y[1], &params);

      params.pinit = y[1];
      // This function is useful if you want to plot, e.g., central
      // pressure vs. total mass. You can also hang onto the run of
      // pressure with radius, which can be interesting when compared to
      // the Newtonian case.
      status = make_grid (params, r, r1, y);
    }

/*
  params.single_star = 0;
  printf("\nnow for a single star!\n");
  printf("%12s %12s %12s %12s\n", "r", "M(r)", "P(r)", "rho(r)");
  // print P and M vs r for Chandrasekhar-mass star
  params.rho_init = 1.0e+13;

  y[1] = pres(params.rho_init, &params); // rho(0) (roughly) for maximum mass
  y[0] = (4.0/3.0) * M_PI * pow(r, 3.0) * rho(y[1], &params);

  params.pinit = y[1];
  status = make_grid(params, r, r1, y);
*/

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  free (eos_tab_pres);
  free (eos_tab_dens);

  eos_tab_dens = NULL;
  eos_tab_pres = NULL;
  fp = NULL;
  pres = NULL;
  rho = NULL;
  return 0;
}


int
make_grid (struct param params, double r, double r1, double y[])
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

  gsl_odeiv2_system sys = { func, fake_jac, 2, &params };
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new
    (&sys, gsl_odeiv2_step_rk8pd,
     1.0e+0, 1.0e-8, 1.0e-8);

  status = gsl_odeiv2_driver_apply (d, &r, r1, y);

  gsl_odeiv2_driver_free (d);
  return 0;
}
