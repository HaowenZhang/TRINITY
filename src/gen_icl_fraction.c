// Calculate the fraction of stellar masses in intra-cluster light
// vs. in central galaxies as a function of halo mass and redshift.
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
#include "check_syscalls.h"
#include "expcache2.h"

extern int64_t num_outputs;
extern struct timestep *steps;

// Calculate the halo mass that corresponds to
// a certain number density at a certain redshift.
double nd_to_mass(double scale, double nd) 
{
  double mass = 17;
  double last_tot = 0;
  double tot = 0;
  for (; tot < nd; mass -= 0.01) {
    last_tot = tot;
    tot += pow(10, mf_cache(scale, mass))/100.0;
  }
  mass += 0.015;
  return (mass - 0.01*(nd-last_tot)/(tot-last_tot)) ;
}



int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int64_t i, j;
  char buffer[1024];

  if (argc<4+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache extension (mcmc output)\n", argv[0]);
    exit(1);
  }
  // Read the model parameters.
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  // Generate the cache arrays to calculate exp10() quickly.
  gen_exp10cache();
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
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);

  
  for (i=0; i<num_outputs; i++) 
  {
    sprintf(buffer, "plots/icl_frac_%f.dat", steps[i].scale);
    // Only look at halos above a certain number density threshold.
    float max_mass = nd_to_mass(steps[i].scale, 1e-6);
    FILE *output = fopen(buffer, "w");
    double total_sm = 0;
    double total_icl = 0;
    for (j=0; j<M_BINS; j++) 
    {
      if (M_MIN+(j+0.5)*INV_BPDEX > max_mass) break;
      total_sm += steps[i].t[j]*steps[i].sm_avg[j];
      total_icl += steps[i].t[j]*steps[i].sm_icl[j];
      fprintf(output, "%f %f %f\n", M_MIN+(j+0.5)*INV_BPDEX, steps[i].icl_frac[j], steps[i].sm_icl[j]/steps[i].sm_avg[j]);
    }
    fprintf(output, "#Avg. ICL frac: %f\n", total_icl/total_sm);
    fclose(output);
  }
  return 0;
}
