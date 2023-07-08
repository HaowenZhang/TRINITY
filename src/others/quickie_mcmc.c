#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/wait.h>
#include "smf.h"
#include "all_smf.h"
#include "observations.h"
#include "jacobi.h"
#include "mt_rand.h"
#include "calc_sfh.h"
#include "mlist.h"
#include "psf_cache.h"

extern float inv_temperature;
extern int64_t num_points;
extern struct smf_fit initial_smf_fit;

void initial_conditions(void);

int main(int argc, char **argv)
{
  int64_t i;
  float temps[7] = {1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0};
  int64_t num_temps = 7;
  init_mcmc_from_args(argc, argv);
  initial_conditions();
  clear_stats();
  set_identity_portion();
  if (CHI2(initial_smf_fit)>1e4) {
    inv_temperature = 0.01; //100
    fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
    mcmc_smf_fit(2000, 0);
    clear_stats();
    inv_temperature = 0.05; //20
    fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
    mcmc_smf_fit(2000, 0);
    stats_to_step(num_points);
  }


  if (CHI2(initial_smf_fit)>1e3) {
    inv_temperature = 0.1; //10
    fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
    clear_stats();
    mcmc_smf_fit(BURN_IN_LENGTH/10.0, 0);
    stats_to_step(num_points);
    clear_stats();
    inv_temperature = 0.15; 
    fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
    mcmc_smf_fit(BURN_IN_LENGTH/10.0, 0);
  }

  clear_stats();
  set_identity_portion();
  inv_temperature = 1;
  fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
  mcmc_smf_fit(BURN_IN_LENGTH*temps[i]/10.0, 0);
  stats_to_step(num_points);
  mcmc_smf_fit(BURN_IN_LENGTH*temps[i]/10.0, 2);
  return 0;
}

