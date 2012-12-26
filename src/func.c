#include <stdio.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "global.h"
#include "eos_pres.h"
#include "eos_rho.h"
#include "param.h"

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
