#ifndef GLOBAL
#define GLOBAL

#ifdef MAIN_FILE
// interpolation machinery
gsl_spline *spline;
gsl_interp_accel *acc;
// bounds for root-finder to invert EOS
// const double rho_root_min = 1.01e+6, rho_root_max = 4.9e+15; // use for bck.eos
// const double rho_root_min = 7.81, rho_root_max = 3.0e+16; // use for eosC
// const double rho_root_min = 1.1e+6, rho_root_max = 9.5e+15; // use for timmes.eos
const double rho_root_min = 1.1e+6, rho_root_max = 9.9e+14;	// use for helmholtz.eos
/* temporary storage of pressure when calling root-finder to invert
   the EOS */
double tmp_pres;
struct param params;
#else
extern gsl_spline *spline;
extern gsl_interp_accel *acc;
// extern const double rho_root_min = 1.01e+6, rho_root_max = 4.9e+15; // use for bck.eos
// extern const double rho_root_min = 7.81, rho_root_max = 3.0e+16; // use for eosC
// extern const double rho_root_min = 1.1e+6, rho_root_max = 9.5e+15; // use for timmes.eos
extern const double rho_root_min, rho_root_max;	// use for helmholtz.eos
extern double tmp_pres;
extern struct param params;

#endif

#endif
