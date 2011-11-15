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

double *eos_tab_pres = NULL, *eos_tab_dens = NULL;

// interpolation stuff
gsl_interp_accel *acc;
gsl_spline *spline;

// root-finding stuff for inverting EOS
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;
gsl_function F;

double tmp_pres;


struct param // parameter list to be passed to ODEs
{
  // pinit: central pressure
  // Gamma: polytropic index
  double Gamma, pinit;
  // single_star: if 0 (true), print grid along r
  int single_star;
  /* temperature. not part of ODE system being integrated so we have
     to leave it in the tag-along struct */
  double temp;
};


int make_grid (struct param params, double r, double r1, double y[]);


// EOS. input: rho, T. output: P
double eos_pres (double rho, struct param *params)
{
  return (gsl_spline_eval (spline, rho, acc));
  //  return (pow(rho, 1.0/params.Gamma)); // polytrope
}


// "reverse" EOS. input: P, T. output: rho
double eos_rho (double pres, struct param *params) {
  int status;
  double rho_low, rho_hi, root;
  tmp_pres = pres;
  do {
    status = gsl_root_fsolver_iterate (s);
    root = gsl_root_fsolver_root (s);
    rho_low = gsl_root_fsolver_x_lower (s);
    rho_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (rho_low, rho_hi, 0, 0.001);
  } while (status == GSL_CONTINUE);
  return (root);
}

double root_func (double rho, void *params) {
  struct param *myparams = (struct param *) params;
  return (tmp_pres - eos_pres(rho, myparams));
}

// the system of ODEs to be integrated
int func (double r, const double y[], double f[], void *params)
{
  struct param *myparams = (struct param *) params;
  double (*pres) (double, struct param *); // EOS pressure pointer
  double (*rho) (double, struct param *); // EOS density pointer
  pres = &eos_pres;
  rho = &eos_rho;

  // for a single star we can print P and M as functions of r
  if (myparams->single_star == 0)
  {
    printf("%12f %12f %12f\n", r, y[0], y[1]);
  }

  if (y[1] < 0.0) // can't have negative pressure
  {
    printf("%12f %12f %12f\n", r, y[0], myparams->pinit);
    return GSL_EBADFUNC; // this flag tells GSL integrator to quit
  } 

  // mass conservation equation
  f[0] = 4.0*M_PI*pow(r, 2.0)*rho(y[1], myparams);

  // TOV equation
  f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.0))
    * (rho(y[1], myparams)
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
    printf("%12f %12f %12f\n", r, y[0], myparams->pinit);
    return GSL_EBADFUNC;
  }

  return GSL_SUCCESS;
}


int main (void)
{
  int i, status;
  struct param params;
  double r, r1, y[2];
  const double pres_min = 1.0e+8, pres_max = 1.0e+45;
  const double rho_min = 1.0e+13, rho_max = 1.0e+17;
  const int MAX = 1000;
  // pressure/density arrays from tabulated EOS data
  int n_eos_pts;
  FILE *fp;
  double (*pres) (double, struct param *); // EOS pressure pointer
  double tmp;

  F.function = &root_func;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  status = gsl_root_fsolver_set (s, &F, pres_min, pres_max);

  pres = &eos_pres;

  fp = fopen("bck.eos", "r");

  if (fp == NULL) {
    fprintf(stderr, "Can't open file!\n");
    exit(1);
  }

  // read in # of EOS data points
  fscanf(fp, "%i", &n_eos_pts);

   printf("# of EOS data points: %i\n", n_eos_pts);
  eos_tab_dens = (double *) malloc (n_eos_pts*sizeof(double));
  eos_tab_pres = (double *) malloc (n_eos_pts*sizeof(double));

 // read in density and pressure data
  for (i = 0; i < n_eos_pts; i++) {
    fscanf(fp, "%le %le %*le %*le", &eos_tab_dens[i], &eos_tab_pres[i]);
   }
  fclose(fp);

  /*
  for (i = 0; i < n_eos_pts; i++) {
    printf ("eos_dens[%i] = %le\t", i, eos_dens[i]);
    printf ("eos_pres[%i] = %le\n", i, eos_pres[i]);
  }
  */

  // set up interpolation machinery for EOS
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, n_eos_pts);
  gsl_spline_init (spline, eos_tab_dens, eos_tab_pres, n_eos_pts);

  printf ("getting ready to interpolate...\n");
  for (i = 0; i < MAX; i++) {
    tmp = rho_min + (double)i*((rho_max-rho_min)/(double) MAX);
    printf("rho = %le; pres = %le\n", tmp, pres(tmp, &params));
  }

  // SKIP THIS STUFF UNTIL EOS (FORWARDS AND BACKWARDS) WORKS!
  /*
  printf ("%12s %12s %12s\n", "R", "M(r=R)", "P(r=0)");

  r = 1.0e+2; // the integrator will go nuts if we start right at r=0
  r1 = 1.0e+15; // some final 'radius' (much larger than actual radius)
  params.single_star = 1;
  for (i = 0; i < MAX; i++)
  {
    central mass always starts at 0-ish
    y[0] = 1.0e-6;
    y[1] = pres(rho_min, params)
      + (double)i*((pres(rho_max, params)
      - pres(rho_min, params))/(double)MAX);
    params.pinit = y[1];
    // This function is useful if you want to plot, e.g., central
    // pressure vs. total mass. You can also hang onto the run of
    // pressure with radius, which can be interesting when compared to
    // the Newtonian case.
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

  */

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free (eos_tab_pres);
  free (eos_tab_dens);
  gsl_root_fsolver_free (s);

  eos_tab_dens = NULL;
  eos_tab_pres = NULL;
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
