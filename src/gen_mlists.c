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

// Total number of snapshots that are calculated in Trinity.
#define TOTAL_OUTPUTS 180
// The redshift of the very first snapshot
#define Z_START 15.0
// The time interval between consecutive snapshots. ***Deprecated***
// because it is too big for the early Universe.
#define STEP_SPACING 80e6
// The minimum halo number density ***per halo mass bin***, below
// which we ignore the halos.
#define ND_CUTOFF -15
#define INV_DZ_CACHE (500.0)
#define DZ_CACHE (1.0/INV_DZ_CACHE)
#define Z_BINS ((int64_t)Z_START*(int64_t)(INV_DZ_CACHE) + 2)

extern struct timestep *steps;
extern int64_t num_outputs;

// The helper cache for calc_step_at_z(). See below.
int64_t stepcache[Z_BINS] = {0};

// Initialize stepcache. stepcache[i] corresponds to a redshift
// z = DZ_CACHE * i, and its value is the index of the snapshot
// with the biggest scale factor that is still smaller than 
// 1 / (1 + z).
void build_stepcache(void) 
{
  int64_t i=0,j;
  double z=0;
  for (j=num_outputs-1; (i<Z_BINS)&&(j>=0); j--) 
  {
    while (i<Z_BINS && (1.0/(1.0+z) > steps[j].scale)) 
    {
      stepcache[i] = j;
      i++;
      z = DZ_CACHE*i;
    }
  }
  // For those without a corresponding snapshot,
  // just assign the very first snapshot.
  while (i<Z_BINS) stepcache[i++] = 0;
}

// Find the closest snapshot number, step, for a given redshift z,
// such that steps[step].scale < 1 / (1 + z). Also calculates
// the distance ***in scale factor dimension*** of 1 / (1 + z)
// from steps[step].scale, f. This is useful when interpolating
// between two consecutive snapshots.
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
  // Find the biggest snapshot number, s, that has steps[s].scale <= 1 / (1 + z)
  while (steps[s].scale > scale && s > 0) s--;
  while (s < num_outputs-1 && (steps[s+1].scale < scale)) s++;

  // If we found the last snapshot, we got no next snapshot to interpolate with.
  // So f should be set to -1 
  if (s==num_outputs-1) *f=-1;

  // If the redshift z too big (i.e., bigger than Z_START) to be covered by any
  // snapshot, also set f to -1.
  else if (s==0 && steps[s].scale > scale) *f=-1;

  // Otherwise, f is the distance between the 1 / (1 + z) and steps[s].scale,
  // normalized by that between steps[s+1].scale and steps[s].scale.
  else *f = (scale - steps[s].scale) / (steps[s+1].scale-steps[s].scale);
  *step = s;
}

// Calculates the correction factor for satellite fraction that corrects 
// the ***old*** total ND (from Appendix G of Behroozi et al. 2013) to 
// the ***new*** total ND (based on the satellite fraction from the
// Universe Machine, Behroozi et al. 2019),
// i.e., ND_new = ND_old * _sat_frac_recorrect().
double _sat_frac_recorrect(double scale, double mass)
{
  double a = scale;
  double a2 = a * a;
  double a3 = a2 * a;
  double a4 = a3 * a;
  double a5 = a4 * a;

  // Calculate the old satellite-to-central ratio from Behroozi et al. (2013)
  double C_a = exp10(-1.91 + 6.23 * a - 15.07 * a2 + 15.02 * a3 - 5.29 * a4);
  double M_cut = 10.66 + 15.93 * a - 21.39 * a2 + 18.2 * a3 - 8.21 * a4;
  double s_to_c_old = C_a * (M_cut - mass);
  if (s_to_c_old < 0) s_to_c_old = 0;

  // Calculate the new satellite-to-central ratio from Behroozi et al. (2019)
  C_a = -0.07378731 + 1.53568848 * a -5.19960262 * a2 +  8.793913 * a3 - 7.30511184 * a4 + 2.35253939 * a5;
  M_cut = 10.57958835 + 16.16331646 * a -24.80597529 * a2 + 19.47148052 * a3 -5.91706418 * a4;
  double s_to_c_new = C_a * (M_cut - mass);
  if (s_to_c_new < 0) s_to_c_new = 0;

  return (1 + s_to_c_new) / (1 + s_to_c_old);
}

// Initialize a certain timestep given its scale factor.
void _init_timestep(struct timestep *t, float scale) {
#define INT_STEPS 8.0
  int64_t i,j;
  float mass;

  memset(t, 0, sizeof(struct timestep));
  for (i=0; i<M_BINS; i++) 
  {
    // Calculate the median halo mass.
    t->med_hm_at_a[i] = M_MIN + (i+0.5)*INV_BPDEX;
    // Further divide each mass bin into INT_STEP sub-bins.
    for (j=0; j<INT_STEPS; j++) 
    {
      //Mass in units of Msun, ***no h^-1***
      mass = M_MIN+(i+(j+0.5)/(INT_STEPS))*INV_BPDEX;
      // Get the halo number density from the cached halo mass functions.
      t->t[i] += pow(10, mf_cache(scale, mass))*INV_BPDEX/INT_STEPS;

    }

    // Ignore halo mass bins with too low number densities.
    if (t->t[i] < pow(10, ND_CUTOFF)) t->t[i] = 0;

    // t->n[i] are the number densities of halos in a certain mass bin that is 
    // inherited from last snapshot. This will be calculated in init_timesteps().
    t->n[i] = 0;
  }
}

// Calculate the integrated merger rate (i.e., the number of mergers per
// central halo, per unit redshift, per unit log10(mass ratio)), given
// the descendant (or merger remnant) mass, mass, and the lower and upper
// mass limits of the satellite halo, m1 and m2, and the scale factor, scale.
double _int_merger_rate(double mass, double m1, double m2, double scale) {
  double dm = (m2-m1)/INT_STEPS;
  double sum = 0;
  int64_t i;

  // Divide the satellite mass range, (m1, m2), into INT_STEPS sub-bins.
  for (i=0; i<INT_STEPS; i++) 
  {
    double m = m1 + (i+0.5)*dm;
    double ratio = m-mass;
    sum += merger_rate_mpeak(mass, ratio, scale);
  }
  return fabs(sum/INT_STEPS);
}

// Calculate the integrated merger rate (i.e., the number of mergers per
// central halo, per unit redshift, per unit log10(mass ratio)), between
// two halo mass bins, b1 and b2, at a certain scale factor, scale.
// Similar to _int_merger_rate() that does the integral over the dimension
// of satellite mass, this function integrates over the dimension of 
// the descendant mass.
double _int2_merger_rate(int64_t b1, int64_t b2, double scale) {
  double dm = INV_BPDEX/INT_STEPS;
  double m_low = b2*INV_BPDEX + M_MIN;
  double m_high = (b2+1)*INV_BPDEX + M_MIN;
  double sum = 0;
  int64_t i;

  // Divide the descendant mass bin, b1, into INT_STEPS sub-bins.
  for (i=0; i<INT_STEPS; i++) 
  {
    double m = b1*INV_BPDEX + M_MIN + (i+0.5)*dm;
    sum += _int_merger_rate(m, m_low, m_high, scale);
  }
  return fabs(sum/INT_STEPS);
#undef INT_STEPS
}

// Calculate the stellar mass loss due to stellar evolutions
// for the n-th snapshot of ts, i.e., ts[n]. Note that the
// stellar mass loss is a function of both the current scale
// factor (ts[n].scale) and the scale factor at which the stars
// formed (ts[k].scale, k=0,1,2...,n).
void calc_smloss(struct timestep *ts, int n) 
{
  int k;
  double a = ts[n].scale;
  // ts[n].smloss[k] is the ***remaining*** fraction of stars
  // that formed between the (k-1)-th and k-th snapshots
  // at the n-th snapshot.
  ts[n].smloss[0] = calc_sm_loss_int(0,ts[0].scale, a);
  for (k=1; k<=n; k++)
    ts[n].smloss[k] = calc_sm_loss_int(ts[k-1].scale,ts[k].scale, a);
}

// Load the timestep information that are pre-calculated and saved.
int64_t load_steps_from_file(char *filename) 
{
  int64_t nout = -1, success = 0, i;
  FILE *input = fopen(filename, "rb");
  if (!input) return 0;
  fread(&nout, sizeof(int64_t), 1, input);
  if (nout == num_outputs) 
  {
    success = fread(steps, sizeof(struct timestep), num_outputs, input);
    if (success == num_outputs) 
    {
      for (i=0; i<num_outputs; i++) 
      {
      	struct timestep *t = steps+i;

        // t->sm_hist, t->icl_stars, and t->sm_inc_hist should be allocated
        // with space and calculated in calc_sfh.c...
      	t->sm_hist = (double *)calloc((num_outputs)*M_BINS, sizeof(double));
      	t->smloss = (double *)calloc(num_outputs, sizeof(double));
      	t->icl_stars = (double *)calloc(num_outputs*M_BINS, sizeof(double));
        t->sm_inc_hist = (double *)calloc(num_outputs*M_BINS, sizeof(double));
        
        // ... but smloss is fixed and can be calculated right here.
      	calc_smloss(steps, i);
      }
    }
  }
  fclose(input);
  return success;
}

// Save the initialized timestep information to the disk.
void write_steps_to_file(char *filename) 
{
  FILE *output = fopen(filename, "wb");
  if (!output) return;
  fwrite(&num_outputs, sizeof(int64_t), 1, output);
  fwrite(steps, sizeof(struct timestep), num_outputs, output);
  fclose(output);
}

// Initialize all the timesteps.
void init_timesteps(void) 
{
  int64_t i,j,k;
  float scale, years_start, dz;
  float mmp_avail[M_BINS];

  // Set the starting year and the number of snapshots.
  years_start = scale_to_years(1.0/(Z_START+1.0));
  num_outputs = TOTAL_OUTPUTS; 
  // allocate the space for all the timesteps
  steps = check_realloc(steps, sizeof(struct timestep)*(num_outputs),
		     "Allocating timesteps.");
  
// #ifdef __APPLE__
  if (load_steps_from_file("steps_PLaexp_sat_frac.bin")) 
  {
    build_stepcache();
    return;
  }
// #endif /* __APPLE__ */

  // The multiplicative factor between two consecutive snapshots' scale factors. 
  double scale_base = pow(1.0 + Z_START, 1.0 / ((double)TOTAL_OUTPUTS - 1.0));
  for (i=0; i<num_outputs; i++) 
  {
    struct timestep *t = steps + i;
    // tp is the previous snapshot.
    struct timestep *tp = steps + i - 1;
    scale = 1.0 / (1.0 + Z_START) * pow(scale_base, i);
    // Initialize snapshot t.
    _init_timestep(t, scale);

    // Allocate space for arrays/interpolation objects that will
    // be used to calculate star formation histories.
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
    t->icl_exp = 0.3*pow(t->scale, -0.7); //Deprecated.
    // dt is the time interval between t and tp.
    t->dt = (i) ? (scale_to_years(t->scale) - scale_to_years(tp->scale)) : 
      (years_start - scale_to_years(0));
    calc_smloss(steps, i);

    // t->c[i] is the number density of halos in the i-th mass bin
    // that newly emerged (instead of inherited from previous
    // snapshots). This will be calculated further below.
    // t->mmp[i*M_BINS+j] is the number of halos in the
    // i-th mass bin that are the most-massive progenitors (MMPs)
    // of the halos in the j-th mass bin. Also calculated
    // further below.
    if (i==0) 
    {
      for (j=0; j<M_BINS; j++) 
      {
	      t->c[j] = t->t[j];
        t->mmp[j+j*M_BINS] = t->t[j];
      }
      continue;
    }

    dz = fabs(1.0/t->scale - 1.0/tp->scale);

    // mmp_avail is the halo number densities that are available
    // to be MMPs for the halos in the current snapshot, or in
    // other words, MMP reservoir. This will be used to prevent 
    // us from taking too many halos as MMPs from any mass bin.
    memcpy(mmp_avail, tp->t, sizeof(float)*M_BINS);

    // Calculating the number density of halos in the k-th mass bin
    // that are merged into the halos in the j-th halo mass bin.
    for (j=0; j<M_BINS; j++)
      for (k=0; k<M_BINS; k++) 
      {
        double mr = t->t[j]*_int2_merger_rate(j, k, scale)*dz*INV_BPDEX;
      	t->merged[k*M_BINS + j] = mr;
        // Remove these merged halos from the reservoir for MMPs.
      	mmp_avail[k] -= mr;
      }

    // After removing all the merged halos (done above),
    // we now count the most-massive progenitors (MMPs)
    // and the newly emerged halos.
    k = M_BINS-1;
    for (j=M_BINS-1; j>=0; j--) 
    {
      // skip the halo mass bins that don't have any halo left
      // to serve as MMPs.
      while (k>=0 && !mmp_avail[k]) k--;

      // The number of newly emerged halos is the total 
      // number density minus those that are inherited from
      // MMPs.
      if (k < 0) { t->c[j] = t->t[j] - t->n[j]; continue; }

      // Find out if the remaining halos in the MMP reservoir
      // can cover the missing part between the total number
      // density and that from other MMPs.
      double missing = t->t[j]-t->n[j];
      // If so, then make up the difference and update
      // t->n[j] = t->t[j].
      if (mmp_avail[k] >= missing) 
      {
      	t->mmp[k*M_BINS+j] = missing;
      	mmp_avail[k] -= missing;
	      t->n[j] = t->t[j];
      }
      // Otherwise, update t->n[j] as much as we can by adding
      // all the available MMPs, and zero-out mmp_avail[k].
      else 
      {
      	t->n[j]+=mmp_avail[k];
      	t->mmp[k*M_BINS+j]=mmp_avail[k];
      	mmp_avail[k]=0;
      	j++;
      }
    }
  }

  build_stepcache();
  write_steps_to_file("steps_PLaexp_sat_frac.bin");
}

