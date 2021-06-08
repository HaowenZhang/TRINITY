// Calculate the specific halo mass accretion rate given
// an input ***intrinsic*** stellar mass, as a function
// of redshift.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "mah.h"

#define MASS_START 7.5
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

// Calculate the halo mass given the stellar mass, sm,
// and the halo--galaxy connection model, c.
double calc_m_at_sm(double sm, double m, struct smf c) 
{
  double sm2 = calc_sm_at_m(m, c);
  double dm_max = 5*INV_BPDEX;
  double dm = dm_max/10.0;
  if (m<=0) {
    m = M_MIN;
    dm_max = 10;
  }
  int64_t count = 0;
  // The secant method.
  while (fabs(sm2-sm)>0.0001 && count<30) 
  {
    double sm3 = calc_sm_at_m(m+dm,c);
    double f = (sm - sm2)/(sm3-sm2);
    double move_m = f*dm;
    if (fabs(move_m) > dm_max)
      return m;
    dm = move_m / 10.0;
    m+=move_m;
    sm2 = calc_sm_at_m(m,c);
    count++;
  }
  if (count>=25) return -1;
  return m;
}

int main(int argc, char **argv)
{
  float z, m;
  struct smf_fit the_smf;
  int i;
  if (argc < 4) 
  {
    fprintf(stderr, "Usage: %s m mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  // Read in the input stellar mass and model parameters.
  m = atof(argv[1]);
  
  // Read in model parameters
  FILE *param_input = check_fopen(argv[3], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  // Fix some model parameters.
  assert_model(&smf);
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[2]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&the_smf);
  
  // Iterate over z=0-15.
  float hm = m + 1.5;
  for (i=0; i<150; i++) 
  {
    z = 0.001 + i*0.1;
    struct smf smf = smhm_at_z(z, the_smf);
    // Calculate the halo mass given the redshift and stellar mass.
    hm = calc_m_at_sm(m, hm, smf);
    // Calculate the corresponding halo mass at z=0
    float mnow = m_now(1.0/(1.0+z), hm);
    // Calculate the average mass accretion rate given the halo mass at z=0.
    float mar = ma_rate_avg(mnow, 1.0/(1.0+z));
    printf("%f %g\n", z, mar/pow(10, hm));
  }

  return 0;
}
