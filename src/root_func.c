#include "eos_pres.h"
#include "global.h"

double
root_func (double rho, void *params)
{
  struct param *root_func_params = (struct param *) params;
  double (*pres) (double, void *);
  pres = &eos_pres;
  return (tmp_pres - pres (rho, root_func_params));
}
