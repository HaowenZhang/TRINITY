#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "mlist.h"
#include "universe_time.h"
#include "observations.h"
#include "check_syscalls.h"
#include "smloss.h"
#include "merger_rates.h"

#define TOTAL_OUTPUTS 180
#define Z_START 15.0
#define STEP_SPACING 80e6
#define ND_CUTOFF -15
#define INV_DZ_CACHE (500.0)
#define DZ_CACHE (1.0/INV_DZ_CACHE)
#define Z_BINS ((int64_t)Z_START*(int64_t)(INV_DZ_CACHE) + 2)

extern struct timestep *steps;
extern int64_t num_outputs;
int64_t stepcache[Z_BINS] = {0};

void build_stepcache(void) {
  int64_t i=0,j;
  double z=0;
  for (j=num_outputs-1; (i<Z_BINS)&&(j>=0); j--) {
    while (i<Z_BINS && (1.0/(1.0+z) > steps[j].scale)) {
      stepcache[i] = j;
      i++;
      z = DZ_CACHE*i;
    }
  }
  while (i<Z_BINS) stepcache[i++] = 0;
}

void calc_step_at_z(double z, int64_t *step, double *f) {
  int64_t s;
  double scale = 1.0/(1.0+z);
  if (!(z>=0)) {
    *step = stepcache[0];
    *f = -1;
    return;
  }
  if (z>Z_START) {
    *step=0;
    *f = -1;
    return;
  }

  s = stepcache[(int64_t)(z*INV_DZ_CACHE)];
  while (steps[s].scale > scale && s > 0) s--;
  while (s < num_outputs-1 && (steps[s+1].scale < scale)) s++;
  if (s==num_outputs-1) *f=-1;
  else if (s==0 && steps[s].scale > scale) *f=-1;
  else *f = (scale - steps[s].scale) / (steps[s+1].scale-steps[s].scale);
  *step = s;
}

double _sat_frac_recorrect(double scale, double mass)
// The correction factor for satellite fraction that corrects the ***old*** total ND to the ***new*** total ND,
// i.e., ND_new = ND_old * _sat_frac_recorrect().
{
  double a = scale;
  double a2 = a * a;
  double a3 = a2 * a;
  double a4 = a3 * a;
  double a5 = a4 * a;
  double C_a = exp10(-1.91 + 6.23 * a - 15.07 * a2 + 15.02 * a3 - 5.29 * a4);
  double M_cut = 10.66 + 15.93 * a - 21.39 * a2 + 18.2 * a3 - 8.21 * a4;
  double s_to_c_old = C_a * (M_cut - mass);
  if (s_to_c_old < 0) s_to_c_old = 0;

  C_a = -0.07378731 + 1.53568848 * a -5.19960262 * a2 +  8.793913 * a3 - 7.30511184 * a4 + 2.35253939 * a5;
  M_cut = 10.57958835 + 16.16331646 * a -24.80597529 * a2 + 19.47148052 * a3 -5.91706418 * a4;
  double s_to_c_new = C_a * (M_cut - mass);
  if (s_to_c_new < 0) s_to_c_new = 0;

  return (1 + s_to_c_new) / (1 + s_to_c_old);
}

void _init_timestep(struct timestep *t, float scale) {
#define INT_STEPS 8.0
  int64_t i,j;
  float mass;

  memset(t, 0, sizeof(struct timestep));
  for (i=0; i<M_BINS; i++) {
    t->med_hm_at_a[i] = M_MIN + (i+0.5)*INV_BPDEX;
    for (j=0; j<INT_STEPS; j++) {
      //Mass in units of Msun (no h)
      mass = M_MIN+(i+(j+0.5)/(INT_STEPS))*INV_BPDEX; //-log10(0.7);
      t->t[i] += pow(10, mf_cache(scale, mass))*INV_BPDEX/INT_STEPS;
      //t->t[i] += pow(10, mf_cache(scale, mass))*INV_BPDEX/INT_STEPS * _sat_frac_recorrect(scale, mass);
    }
    if (t->t[i] < pow(10, ND_CUTOFF)) t->t[i] = 0;
    t->n[i] = 0;
  }
}


double _int_merger_rate(double mass, double m1, double m2, double scale) {
  double dm = (m2-m1)/INT_STEPS;
  double sum = 0;
  int64_t i;
  for (i=0; i<INT_STEPS; i++) {
    double m = m1 + (i+0.5)*dm;
    double ratio = m-mass;
    sum += merger_rate_mpeak(mass, ratio, scale);
  }
  return fabs(sum/INT_STEPS);
}

double _int2_merger_rate(int64_t b1, int64_t b2, double scale) {
  double dm = INV_BPDEX/INT_STEPS;
  double m_low = b2*INV_BPDEX + M_MIN;
  double m_high = (b2+1)*INV_BPDEX + M_MIN;
  double sum = 0;
  int64_t i;
  for (i=0; i<INT_STEPS; i++) {
    double m = b1*INV_BPDEX + M_MIN + (i+0.5)*dm;
    sum += _int_merger_rate(m, m_low, m_high, scale);
  }
  return fabs(sum/INT_STEPS);
#undef INT_STEPS
}

void calc_smloss(struct timestep *ts, int n) {
  int k;
  double a = ts[n].scale;
  ts[n].smloss[0] = calc_sm_loss_int(0,ts[0].scale, a);
  for (k=1; k<=n; k++)
    ts[n].smloss[k] = calc_sm_loss_int(ts[k-1].scale,ts[k].scale, a);
}


int64_t load_steps_from_file(char *filename) {
  int64_t nout = -1, success = 0, i;
  FILE *input = fopen(filename, "rb");
  if (!input) return 0;
  fread(&nout, sizeof(int64_t), 1, input);
  if (nout == num_outputs) {
    success = fread(steps, sizeof(struct timestep), num_outputs, input);
    if (success == num_outputs) {
      for (i=0; i<num_outputs; i++) {
	struct timestep *t = steps+i;
	t->sm_hist = (double *)calloc((num_outputs)*M_BINS, sizeof(double));
	t->smloss = (double *)calloc(num_outputs, sizeof(double));
	t->icl_stars = (double *)calloc(num_outputs*M_BINS, sizeof(double));
  t->sm_inc_hist = (double *)calloc(num_outputs*M_BINS, sizeof(double));
	// t->spline = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
	// t->spline_sfr = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
	calc_smloss(steps, i);
      }
    }
  }
  fclose(input);
  return success;
}

void write_steps_to_file(char *filename) {
  FILE *output = fopen(filename, "wb");
  if (!output) return;
  fwrite(&num_outputs, sizeof(int64_t), 1, output);
  fwrite(steps, sizeof(struct timestep), num_outputs, output);
  fclose(output);
}

void init_timesteps(void) {
  int64_t i,j,k;
  float scale, years_start, dz;
  float mmp_avail[M_BINS];

  years_start = scale_to_years(1.0/(Z_START+1.0));
  num_outputs = TOTAL_OUTPUTS; // num_outputs = ceil(-years_start/STEP_SPACING)+1;
  steps = check_realloc(steps, sizeof(struct timestep)*(num_outputs),
		     "Allocating timesteps.");
  
#ifdef __APPLE__
  if (load_steps_from_file("steps_PLaexp_sat_frac.bin")) {
    build_stepcache();
    return;
  }
#endif /* __APPLE__ */

  double scale_base = pow(1.0 + Z_START, 1.0 / ((double)TOTAL_OUTPUTS - 1.0));
  for (i=0; i<num_outputs; i++) {
    struct timestep *t = steps + i;
    struct timestep *tp = steps + i - 1;
    scale = 1.0 / (1.0 + Z_START) * pow(scale_base, i); // scale = years_to_scale(years_start + i*STEP_SPACING);
    _init_timestep(t, scale);
    t->sm_hist = (double *)calloc((num_outputs)*M_BINS, sizeof(double));
    t->smloss = (double *)calloc(num_outputs, sizeof(double));
    t->icl_stars = (double *)calloc(num_outputs*M_BINS, sizeof(double));
    t->sm_inc_hist = (double *)calloc(num_outputs*M_BINS, sizeof(double));
    t->spline = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
    t->spline_sfr = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
    t->spline2 = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
    t->spline_sfr2 = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
    // t->spline_uv = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
    // t->spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, M_BINS);
    t->scale = scale;
    t->icl_exp = 0.3*pow(t->scale, -0.7);
    t->dt = (i) ? (scale_to_years(t->scale) - scale_to_years(tp->scale)) : 
      (years_start - scale_to_years(0)); // t->dt = (i) ? STEP_SPACING : (years_start - scale_to_years(0)); 
    calc_smloss(steps, i);

    if (i==0) {
      for (j=0; j<M_BINS; j++) {
	t->c[j] = t->t[j];
        t->mmp[j+j*M_BINS] = t->t[j];
      }
      continue;
    }

    dz = fabs(1.0/t->scale - 1.0/tp->scale);
    memcpy(mmp_avail, tp->t, sizeof(float)*M_BINS);
    for (j=0; j<M_BINS; j++)
      for (k=0; k<M_BINS; k++) {
	double mr = t->t[j]*_int2_merger_rate(j, k, scale)*dz*INV_BPDEX;
	t->merged[k*M_BINS + j] = mr;
	mmp_avail[k] -= mr;
      }

    k = M_BINS-1;
    for (j=M_BINS-1; j>=0; j--) {
      while (k>=0 && !mmp_avail[k]) k--;
      if (k < 0) { t->c[j] = t->t[j] - t->n[j]; continue; }
      double missing = t->t[j]-t->n[j];
      if (mmp_avail[k] >= missing) {
	t->mmp[k*M_BINS+j] = missing;
	mmp_avail[k] -= missing;
	t->n[j] = t->t[j];
      }
      else {
	t->n[j]+=mmp_avail[k];
	t->mmp[k*M_BINS+j]=mmp_avail[k];
	mmp_avail[k]=0;
	j++;
      }
    }
  }

  // // Apply the correction on the merger rates (decrease by 50%) 
  // // and the halo mass function (increase the satellites by 25%)
  // for (i=0; i<num_outputs; i++)
  // {
  //   struct timestep *t = steps + i;
  //   struct timestep *tp = steps + i - 1;
  //   double a = t->scale;

  //   double C_a = exp10(-1.91 + 6.23 * a - 15.07 * a * a + 15.02 * a * a * a - 5.29 * a * a * a * a);
  //   double M_cut = 10.66 + 15.93 * a - 21.39 * a * a + 18.2 * a * a * a - 8.21 * a * a * a * a;
  //   double ap, C_ap, M_cutp;
    
  //   if (i==0)
  //   {
  //     C_ap = C_a;
  //     M_cutp = M_cut;    
  //   }

  //   else
  //   {
  //     ap = tp->scale;
  //     C_ap = exp10(-1.91 + 6.23 * ap - 15.07 * ap * ap + 15.02 * ap * ap * ap - 5.29 * ap * ap * ap * ap);
  //     M_cutp = 10.66 + 15.93 * ap - 21.39 * ap * ap + 18.2 * ap * ap * ap - 8.21 * ap * ap * ap * ap;
  //   }


  //   for (j=0; j<M_BINS; j++)
  //   {
  //     double m = M_MIN + (j + 0.5) * INV_BPDEX;
  //     double s_to_c = C_a * (m - M_cut); //The number ratio of satellite and central galaxies.
  //     // add 25% more satellite galaxies to each mass bin, since they haven't been merged yet.
  //     t->t[j] = t->t[j] / (1 + s_to_c) + t->t[j] * s_to_c / (1 + s_to_c) * (1 + 0.25);
  //     t->n[j] = t->n[j] / (1 + s_to_c) + t->n[j] * s_to_c / (1 + s_to_c) * (1 + 0.25);
  //     t->c[j] = t->c[j] / (1 + s_to_c) + t->c[j] * s_to_c / (1 + s_to_c) * (1 + 0.25);

  //     for (k=0; k<M_BINS; k++)
  //     {
  //       t->merged[k*M_BINS + j] *= 0.5;
  //       double mk = M_MIN + (k + 0.5) * INV_BPDEX; //The halo mass for the k-th bin.
  //       double s_to_c_k = C_ap * (mk - M_cutp); //The satellite to central number ratio. Since we're dealing with MMP, we should use the fraction for the previous snapshot.
  //       t->mmp[k*M_BINS + j] = t->mmp[k*M_BINS + j] / (1 + s_to_c_k) + t->mmp[k*M_BINS + j] * s_to_c_k / (1 + s_to_c_k) * (1 + 0.25);
  //     }
  //   }
  // }
  
  build_stepcache();
  write_steps_to_file("steps_PLaexp_sat_frac.bin");
}

