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
extern float z_max, z_min;

void initial_conditions(void);

int main(int argc, char **argv)
{
  int64_t i;
  double chi2[2000];
  float temps[7] = {1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0};
  int64_t num_temps = 7;
  init_mcmc_from_args(argc, argv);
  initial_conditions();
  random_step(&initial_smf_fit);
  float delta_0 = DELTA(initial_smf_fit);
  for (i=0; i<2000; i++) {
    DELTA(initial_smf_fit) = delta_0 + delta_0*(i-1000)/1000.0;
    chi2[i] = all_smf_chi2_err(initial_smf_fit);
  }
  double total = 0;
  double best_delta = delta_0 - 1.0;
  int64_t best_i = 0;
  for (i=0; i<2000; i++) {
    total += exp(-0.5*chi2[i]);
    if (chi2[best_i] > chi2[i]) { best_i = i; best_delta = delta_0 + delta_0*(i-1000)/1000.0; }
  }
  double sum=0;
  //printf("Total: %e\n", total);
  //for (i=0; i<2000; i++) chi2[i] += log(total)/2.0;
  for (i=0; i<2000; i++) {
    sum += exp(-0.5*chi2[i]);
    if (sum/total > 0.685) break;
  }
  //  printf("Initial i: %"PRId64"\n", i);
  int64_t i_start = 0;
  int64_t i_stop = i;
  int64_t best_i_start = 0;
  int64_t best_i_stop = i;
  for (i_start = 1; i_stop < 2000 && i_start < 1999; i_start++) {
    sum -= exp(-0.5*chi2[i_start]);
    while (sum/total < 0.685 && i_stop < 1999) {
      i_stop++;
      sum += exp(-0.5*chi2[i_stop]);
    }
    if (sum/total < 0.685) break;
    if ((i_stop-i_start) < (best_i_stop-best_i_start)) {
      best_i_stop = i_stop;
      best_i_start = i_start;
    }
  }
  printf("%f %e %e %e %e\n", 0.5*(z_max+z_min), best_delta, (best_i_stop-best_i)*delta_0/1000.0, (best_i-best_i_start)*delta_0/1000.0, chi2[best_i]);
  return 0;
}

