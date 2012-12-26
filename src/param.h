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
