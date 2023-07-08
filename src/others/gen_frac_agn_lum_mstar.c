// Generate quasar luminosity functions broken up into
// the contributions from different galaxy mass bins,
// as functions of redshift.
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

#define MSTAR_MIN 8 //minimum/maximum stellar mass for intrinsic stellar mass functions.
#define MSTAR_MAX 12
#define MSTAR_BPDEX 10 // # of bins per dex in stellar mass
#define MSTAR_BINS ((int64_t)((MSTAR_MAX-MSTAR_MIN)*MSTAR_BPDEX + 2))
#define MSTAR_INV_BPDEX (1.0/((double)MSTAR_BPDEX))

int main(int argc, char **argv)
{
  int64_t i, j, k;
  struct smf_fit smf;
  double l, mstar;
  if (argc<4+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s z lbol_min mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);
  double lbol_min = atof(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+4]);

  nonlinear_luminosity = 1;
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[3]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  int64_t step;
  double f;
  // Calculate the # of snapshot that is the closest to z.
  calc_step_at_z(z, &step, &f);
  double inv_scatter_sm = 1.0 / steps[step].smhm.scatter;
  double inv_scatter_mbh = 1.0 / steps[step].smhm.bh_scatter;
  double gauss_norm_sm = 1 / sqrt(2 * M_PI) * inv_scatter_sm;
  double gauss_norm_mbh = 1 / sqrt(2 * M_PI) * inv_scatter_mbh;

  // Culumative fraction of BHs above different Eddington ratios,
  // in each halo mass bin. This is done so that we don't have to 
  // integrate everytime in the following calculations. Note that
  // these fractions are among ***active BHs*** only, so an additional
  // duty cycle factor (steps[step].bh_duty[j]) is needed when using
  // these pre-calculated fractions.
  double bher_cdf[M_BINS*(BHER_BINS+1)] = {0};
  for (i=0; i<M_BINS; i++)
  {
    // bher_cdf[i*(BHER_BINS+1)+BHER_BINS] = steps[step].bher_dist[i*BHER_BINS+BHER_BINS-1];
    for (j=BHER_BINS-1; j>-1; j--)
      bher_cdf[i*(BHER_BINS+1)+j] = steps[step].bher_dist[i*BHER_BINS+j] + bher_cdf[i*(BHER_BINS+1)+j+1];
    for (j=BHER_BINS-1; j>-1; j--)
      bher_cdf[i*(BHER_BINS+1)+j] /= steps[step].bher_dist_norm[i];
  }


  // Calculate the fractions of galaxies hosting active BHs (above the input bolometric luminosity limit)
  // as a function of galaxy mass.
  for (i=0; i<MSTAR_BINS; i++) 
  {
    
    // Calculate the median BH mass for this galaxy mass.
    double sm = MSTAR_MIN + i * MSTAR_INV_BPDEX;
    // double bm = bulge_mass(sm + steps[step].smhm.mu, steps[step].scale);
    double bm = bulge_mass(sm, steps[step].scale);
    double mbh_med = calc_bh_at_bm(bm, steps[step].smhm);

    // number densities of ALL SMBHs and SMBHs above the luminosity limit.
    double nd_tot = 0;
    double nd_above_lum = 0;


    for (j=0; j<M_BINS; j++)
    {
      if (!steps[step].t[j]) continue;

      // double arg = exp(-1.12882 * (steps[step].log_sm[j] + steps[step].smhm.mu - 10.1993));
      // double k_scatter = 1  + 1.12882 / M_LN10 * arg / (1 + arg);

      // double inv_scatter_sm = 1.0 / steps[step].smhm.scatter * k_scatter;
      // double gauss_norm_sm = 1 / sqrt(2 * M_PI) * inv_scatter_sm;



      // Calculate the probability of this halo mass bin's hosting this galaxy mass
      double dsm = (sm - steps[step].smhm.mu - steps[step].log_sm[j]) * inv_scatter_sm;
      double prob_sm = gauss_norm_sm * exp(-0.5 * dsm * dsm) * MSTAR_INV_BPDEX;
      double log_bh_focc = log10(steps[step].bh_f_occ[j]);

      for (k=0; k<MBH_BINS; k++)
      {
        // Calculate the probability of this BH mass hosted by this galaxy
        double mbh = MBH_MIN + (k + 0.5) * MBH_INV_BPDEX;
        double dmbh = (mbh - mbh_med + log_bh_focc) * inv_scatter_mbh;
        double prob_mbh = gauss_norm_mbh * exp(-0.5 * dmbh * dmbh) * MBH_INV_BPDEX;

        // Calculate the Eddington ratio limit corresponding to the luminosity limit.
        double bh_eta = lbol_min - 38.1 - mbh;
        double eta_frac = bh_eta - (steps[step].bh_eta[j] + log10(steps[step].bh_f_occ[j]) 
          + (steps[step].smhm.rho_bh - 1) * (mbh - steps[step].log_bh_mass[j] - steps[step].smhm.log_bh_scatter_corr));

        double f_eta = (eta_frac - steps[step].ledd_min[j]) * steps[step].ledd_bpdex[j];
        int64_t b_eta = f_eta; f_eta -= b_eta;
        if (b_eta < 0) {b_eta = f_eta = 0;}
        else if (b_eta >= BHER_BINS) {b_eta = BHER_BINS-1; f_eta=1;}

        // Get the cumulative fraction from the pre-calculated CDFs
        double frac_above_lum = bher_cdf[j*(BHER_BINS+1)+b_eta] + 
                        f_eta * (bher_cdf[j*(BHER_BINS+1)+b_eta+1] - bher_cdf[j*(BHER_BINS+1)+b_eta]);

        // Update the number densities of all SMBHs(galaxies) and active SMBHs.
        nd_tot += prob_sm * prob_mbh * steps[step].t[j] * steps[step].bh_f_occ[j];
        nd_above_lum += prob_sm * prob_mbh * steps[step].t[j] * steps[step].bh_duty[j] * frac_above_lum;

      }

    }

    printf("%f %e %e %e\n", sm, nd_tot, nd_above_lum, nd_above_lum/nd_tot); //, log10(calc_quasar_lf(l+6.25, z))-2.5);
    
    
  }

  return 0;
}

