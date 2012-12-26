#include <gsl/gsl_odeiv2.h>
#include "func.h"
#include "make_grid.h"

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
