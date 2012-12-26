#define MAIN_FILE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include "global.h"
#include "eos_pres.h"
#include "eos_rho.h"
#include "param.h"

int
main (void)
{
  spline = NULL;
  acc = NULL;
  int i, status;
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
  fp = fopen ("../EOS/helmholtz.eos", "r");

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
