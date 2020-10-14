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

#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }


int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);

  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  double t,m;
  printf("#1+z dt M_h M_bh M_bh_old Mbh_merge_m Mbh_merge_w Mbh_unmerged sm Mbh_old_self Mbh_old_next\n");
  for (t=0; t<num_outputs-1; t+=1.0/3.0) {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double bh_merge_mchar = (1.0 - f) * steps[i].smhm.bh_merge + f * steps[i+1].smhm.bh_merge;
     double bh_merge_width = (1.0 - f) * steps[i].smhm.bh_merge_width + f * steps[i+1].smhm.bh_merge_width;
     double dt = (1.0 - f) * steps[i].dt + f * steps[i+1].dt;
    for (m=8; m<15; m+=0.05) {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double bh_mass, bh_mass_old, bh_unmerged, log_bh_acc_rate, log_sm, bh_mass_old_self, bh_mass_old_next;
      LINTERP(bh_mass,bh_mass_avg);
      LINTERP(bh_mass_old, old_bh_mass);
      LINTERP(bh_mass_old_next, old_bh_mass_next);
      LINTERP(bh_mass_old_self, old_bh_mass_self);
      LINTERP(log_bh_acc_rate,bh_acc_rate);
      LINTERP(log_sm, sm_avg);
      // INTERP(bh_merge_rate,bh_merge_rate);
      LINTERP(bh_unmerged,bh_unmerged);
      if (!isfinite(log_bh_acc_rate)) continue;
      
      // bh_merge_rate = (bh_merge_rate > 1e-5) ? log10(bh_merge_rate) : -5;
      // bh_unmerged = (bh_unmerged > 0.1 ) ? log10(bh_unmerged) : -1;
      // bh_unmerged = log10(bh_unmerged);
  // float years = scale_to_years(1.0 / zp1);
      printf("%f %f %f %f %f %f %f %f %f %f %f\n", zp1, dt, m, bh_mass, bh_mass_old, bh_merge_mchar, bh_merge_width, bh_unmerged, log_sm, bh_mass_old_self, bh_mass_old_next);
    }
  }
  return 0;
}
