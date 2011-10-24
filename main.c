/*

Copyright 2011 Brian Friesen. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other materials
   provided with the distribution.

THIS SOFTWARE IS PROVIDED BY BRIAN FRIESEN ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Brian Friesen.

*/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <math.h>

struct param // parameter list to be passed to ODEs
{
  // pinit: central pressure
  // Gamma: polytropic index
  double Gamma, pinit;
};

/* EOS function so we don't have to hard-code it into the structure of the
   ODEs. For now, since the 2 ODEs are for P and M, we'll assume the EOS is
   in the form rho = rho(P), which, I think, is usually the inverse of the
   structure of most real EOS routines. I'm doing it this way just so I
   don't have to change one of the ODE variables from P to rho. In the long
   run that would be the better way to do this, just so that it would play
   nicer with external EOSs. */
double eos_rho(double pres, double Gamma)
{
  return (pow(pres, 1.0/Gamma)); // polytrope
}

int func (double r, const double y[], double f[], void *params)
{
  struct param myparams = *(struct param *)params;
  double (*rho)(double, double); // EOS pointer
  rho = &eos_rho;

  if (y[1] < 0.0) // can't have negative pressure
  {
    printf("%12f %12f %12f\n", r, y[0], myparams.pinit);
    return GSL_EBADFUNC; // this flag tells GSL integrator to quit
  } 

  // mass conservation equation
  f[0] = 4.0*M_PI*pow(r, 2.0)*rho(y[1], myparams.Gamma);

  // TOV equation
  f[1] = -(1.0/pow(r, 2.0)) * (rho(y[1], myparams.Gamma) + y[1]) *
         (y[0] + 4.0*M_PI*pow(r, 3.0)*y[1]) / (1.0 - (2.0*y[0]/r));

// ODEs tend to return NaNs at the same radius where P < 0
  if ((gsl_isnan(f[0]) || gsl_isnan(f[1]))) 
  {
    /* print values at the radius where ODEs diverge. this is almost always
       at the surface */
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

  params.Gamma = 5.0/3.0; // electron-degenerate EOS
//  params.Gamma = 3.0; // neutron-degenerate EOS

  printf("%12s %12s %12s\n", "R", "M(r=R)", "P(r=0)");
  r = 1.0e-5; // the integrator will go nuts if we start right at r=0
  r1 = 3.0; // some final 'radius' (much larger than actual radius)
  for (i = 0; i < MAX; i++)
  {
    // central mass always starts at 0 (or very close, for numerical reasons)
    y[0] = 1.0e-6;
    y[1] = pmin + (double)i*((pmax - pmin)/(double)MAX);
    params.pinit = y[1];
    /* This function is useful if you want to plot, e.g., central
     * pressure vs. total mass. You can also hang onto the run of pressure
     * with radius, which can be interesting when compared to the Newtonian
     * case. */
    status = make_grid(params, r, r1, y);
  }
  
  return 0;
}

int make_grid (struct param params, double r, double r1, double y[])
{
  int status;
  int *fake_jac;
 /* The integrator RK8PD (Runge-Kutta Prince-Dormand) doesn't need the jacobian
  * (some more sophisticated integrators do) so make this a null pointer.
  * The linker will complain that the pointer type is incompatable since it
  * doesn't point to a function with a bunch of arguments, but the pointer
  * won't be used anyway so it doesn't matter. */
  fake_jac = 0;
		   
  gsl_odeiv2_system sys = {func, fake_jac, 2, &params};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new
                         (&sys, gsl_odeiv2_step_rk4, 1.0e-6, 1.0e-6, 0.0);

  status = gsl_odeiv2_driver_apply(d, &r, r1, y);

  gsl_odeiv2_driver_free (d);
  return 0;
}
