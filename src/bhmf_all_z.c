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
#include "expcache2.h"

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  double m;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  // double z = atof(argv[1]);
  double z;
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);

  double zs[12] = {0};
  for (i=2; i<12; i++) zs[i] = i - 1;
  zs[0] = 0.1; zs[1] = 0.5;

  fprintf(stdout, "#z BH_Mass ND ND(Active)\n");
  for (i=0; i<12; i++)
  {
    int64_t step;
    double f;
    calc_step_at_z(zs[i], &step, &f);
    double s1 = steps[step].smhm.scatter;
    double s2 = steps[step].smhm.bh_scatter;
    double s = sqrt(s1*s1+s2*s2);
    for (m=4; m<10.5; m+=0.1) 
    {
      fprintf(stdout, "%f %f %e %e\n", zs[i], m, calc_bhmf(m, zs[i]), calc_active_bhmf(m,zs[i]));
    }
  }

  //// Redo the calculation, but now without merger contributions.
  remove_bh_merger();
  fprintf(stderr, "#z BH_Mass ND ND(Active)\n");
  for (i=0; i<12; i++)
  {
    int64_t step;
    double f;
    calc_step_at_z(zs[i], &step, &f);
    double s1 = steps[step].smhm.scatter;
    double s2 = steps[step].smhm.bh_scatter;
    double s = sqrt(s1*s1+s2*s2);
    for (m=4; m<10.5; m+=0.1) 
    {
      fprintf(stderr, "%f %f %e %e\n", zs[i], m, calc_bhmf(m, zs[i]), calc_active_bhmf(m,zs[i]));
    }
  }

  
  
  // printf("#double power-law norm: %e\n", doublePL_norm(-0.6, -10, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -10, 2, NULL));
  // printf("#double power-law norm: %e\n", doublePL_norm(-0.6, -4, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -4, 2, NULL));
  // printf("#BH_Alpha: %f\n", steps[step].smhm.bh_alpha);
  // printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  
  // printf("#BH_Scatter at fixed mass: %f\n", s);
  
  // FILE * pfile;
  // pfile = fopen('./ABHMF_check.txt', 'w');
  // fprintf(pfile, "Mbh Mh ND f_active\n");
  
  // fclose(pfile);
  return 0;
}
