#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <sys/stat.h>
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
#include "stringparse.h"
#include "mt_rand.h"
// #include "inthash.h"

extern int64_t num_outputs;
extern struct timestep *steps;

#define A_UV(x) ((x).t_tdyn)
struct catalog_halo {
  int64_t id, descid, upid;
  int32_t flags; float uparent_dist;
  float pos[6], vmp, lvmp, mp, m, v, r;
  float rank1, rank2, ra, rarank, t_tdyn;
  float sm, icl, sfr, obs_sm, obs_sfr, obs_uv;
};

#define IGNORE_FLAG 16
// The indices to locate scale factors in the file name strings.
// These indices only apply to files like "sfr_catalog_1.002310.bin".
#define SCALE_START -12
#define SCALE_LENGTH 8

struct catalog_halo *read_from_um_catalog(char *filename, int64_t *num_halos)
{
  struct stat file_stat;
  FILE *input=NULL;
  struct catalog_halo *halos = NULL;
  int64_t i;

  input = fopen(filename, "rb");
  if (!input || fstat(fileno(input), &file_stat)<0) {
    fprintf(stderr, "[Error] Could not open/stat file %s!\n", filename);
    exit(EXIT_FAILURE);
  }

  if (file_stat.st_size % sizeof(struct catalog_halo)) {
    fprintf(stderr, "[Error] Size of file %s not evenly divisible by halo structure size (%"PRId64")\n",
      filename, (int64_t)sizeof(struct catalog_halo));
    exit(EXIT_FAILURE);
  }
  
  *num_halos = file_stat.st_size /  sizeof(struct catalog_halo);
  halos = malloc(sizeof(struct catalog_halo)*(*num_halos));
  if (!halos) {
    fprintf(stderr, "[Error] Could not allocate memory for halos.");
    exit(EXIT_FAILURE);
  }
  
  if (fread(halos, sizeof(struct catalog_halo), *num_halos, input)<0) {
    fprintf(stderr, "[Error] Could not read file %s!\n", filename);
    exit(EXIT_FAILURE);    
  }
  fclose(input);

  return halos;
}


float biterp (float a, float b, float c, float d, float f1, float f2) {
  float al = log10(a);
  float bl = log10(b);
  float cl = log10(c);
  float dl = log10(d);
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return (e+f2*(f-e));
}


void add_agn_to_catalog(struct smf_fit f, char *filename) {
  char buffer[1024];
  FILE *input, *output;
  int n, ts, i;
  float sm, obs_sm;
  float mbh; 
  double p_eta[4] = {0};
  double p_lbol[5] = {0};
  double eta_aird = -1; //The selection criteria from Aird2019.
  double lbol_aird = 45;
  double p_aird = 0;
  char scale_string[1024];
  int64_t filename_len = strlen(filename);
  strncpy(scale_string, filename + filename_len + SCALE_START, SCALE_LENGTH);

  double scale = atof(scale_string);
  fprintf(stderr, "scale_string=%s, scale=%f\n", scale_string, scale);
  

  int64_t num_halos = 0;
  struct catalog_halo *halos = read_from_um_catalog(filename, &num_halos);
  fprintf(stderr, "num_halos=%d\n", num_halos);
  
  strncpy(buffer, "Catalogs/final_eff_z/sm_agn_catalog_a", 1024);
  strcat(buffer, scale_string);
  strcat(buffer, ".dat");
  fprintf(stderr, "Output file name: %s\n", buffer);  
  n = strlen(buffer);


  output = check_fopen(buffer, "w");
  fprintf(output, "#scale(0) id(1) desc_id(2) upid(3) flags(4) uparent_dist(5) x(6) y(7) z(8) vx(9) vy(10) vz(11) mvir(12) vmax(13) Mpeak(14) Vpeak(15) rvir(16) Rank1(17) Rank2(18) A_UV(19) Log(SM)(20) Log(ICL)(21) log(SFR)(22) Log(obs_SM)(23) Log(obs_SFR)(24) obs_UV(25) Mbh(26) p_eta_0001(27) p_eta_001(28) p_eta_01(29) p_eta_1(30) p_lbol_41(31) p_lbol_42(32) p_lbol_43(33) p_lbol_44(34) p_lbol_45(35) p_eta01_lbol45(36)\n");

  
  // Header lines
  fprintf(output, "#Field explanations:\n");
  fprintf(output, "#**Note that halo masses are in Msun/h and stellar masses/SFRs are in Msun.\n");
  fprintf(output, "#Scale: scale factor\n");
  fprintf(output, "#ID: Unique halo ID\n");
  fprintf(output, "#DescID: ID of descendant halo (or -1 at z=0).\n");
  fprintf(output, "#UPID: -1 for central halos, otherwise, ID of largest parent halo\n");
  fprintf(output, "#Flags: Mostly internal UniverseMachine info.  However, halos with bit 4 set in the flags (i.e., flags & (2**4) is true) should be ignored.\n");
  fprintf(output, "#Uparent_Dist: Ignore\n");
  fprintf(output, "#X Y Z: halo position (comoving Mpc/h)\n");
  fprintf(output, "#VX VY VZ: halo velocity (physical peculiar km/s)\n");
  fprintf(output, "#M: Halo mass (Bryan & Norman 1998 virial mass, Msun/h)\n");
  fprintf(output, "#V: Halo vmax (physical km/s)\n");
  fprintf(output, "#MP: Halo peak historical mass (BN98 vir, Msun/h)\n");
  fprintf(output, "#VMP: Halo vmax at the time when peak mass was reached.\n");
  fprintf(output, "#R: Halo radius (BN98 vir, comoving kpc/h)\n");
  fprintf(output, "#Rank1: halo rank in Delta_vmax (see UniverseMachine paper)\n");
  fprintf(output, "#Rank2, RA, RARank: Ignore\n");
  fprintf(output, "#A_UV: UV attenuation (mag)\n");
  fprintf(output, "#SM: True stellar mass (Msun)\n");
  fprintf(output, "#ICL: True intracluster stellar mass (Msun))\n");
  fprintf(output, "#SFR: True star formation rate (Msun/yr)\n");
  fprintf(output, "#Obs_SM: observed stellar mass, including random & systematic errors (Msun)\n");
  fprintf(output, "#Obs_SFR: observed SFR, including random & systematic errors (Msun/yr)\n");
  fprintf(output, "#Obs_UV: Observed UV Magnitude (M_1500 AB)\n");

  fprintf(output, "#Mbh: SMBH mass (Msun)\n");
  fprintf(output, "#P_eta_${eta_value}: The cumulative probability that the galaxy hosts an AGN with an Eddington ratio above eta_value.\n");
  fprintf(output, "#P_lbol_${log(lbol)_value}: The cumulative probability that the galaxy hosts an AGN with a bolometric luminosity above log(lbol)_value.\n");
  fprintf(output, "#P_eta${eta_value}_lbol${log(lbol)_value}: The cumulative probability that the galaxy hosts an AGN with an Eddington ratio above eta_value, AND a bolometric luminosity above log(lbol)_value.\n");

  ts = num_outputs-2;

  for (i=0; i<(num_halos); i++)
  //for (i=0; i<(1000); i++)
  {
    // Ignore the halos that are already ignored in the Universe Machine catalogs.
    if (halos[i].flags & IGNORE_FLAG) continue;

    // Find the snapshot that is the closest to the Universe Machine catalog in the scale factor space.
    while (steps[ts].scale > scale && ts>0) ts--;
    if (!ts) break;

    // Determine how much we should interpolate the model parameters.
    float tf = (scale - steps[ts].scale) / (steps[ts+1].scale - steps[ts].scale);
    
    // Ignore the halos with zero halo mass.
    if (!(halos[i].mp>0) || !(halos[i].m>0)) continue;
    float logm = log10(halos[i].mp/0.678);

    // Calculate stellar mass based on Trinity.
    float sm1 = calc_sm_at_m(logm, steps[ts]);
    float sm2 = calc_sm_at_m(logm, steps[ts+1]);
    sm = sm1 + tf*(sm2-sm1);
    
#define INTERP(y) (steps[ts].smhm.y + tf*(steps[ts+1].smhm.y - steps[ts].smhm.y))
    // Interpolate to get relevant model parameters
    float ins_scatter = INTERP(scatter);
    float obs_scatter = INTERP(obs_scatter);
    float rho = INTERP(sfr_sm_corr);
    float mu = INTERP(mu);
    float kappa = INTERP(kappa);

    float bh_gamma = INTERP(bh_gamma);
    float bh_beta = INTERP(bh_beta);
    float bh_scatter = INTERP(bh_scatter);
    float bh_duty = INTERP(bh_duty);
    float dc_mbh = INTERP(dc_mbh);
    float dc_mbh_w = INTERP(dc_mbh_w);

    float sm_offset = normal_random(0, 1);
    float sfr_offset = normal_random(0, 1);
    // The random scatter in the Mbh--Mbulge relation.
    float mbh_offset = normal_random(0, 1);

    sm = sm + sm_offset*ins_scatter;
    obs_sm = sm + normal_random(0, obs_scatter) + mu;

    // Calculate the BH mass
    // float bm = bulge_mass(log10(halos[i].obs_sm), scale);
    float bm = bulge_mass(obs_sm, scale);
    mbh = bh_beta + bh_gamma * (bm - 11.0);
    mbh += mbh_offset * bh_scatter;

    // Calculate the duty cycle.
    double f_mass = exp((mbh - dc_mbh) / dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    // f_mass = f_mass < 1? f_mass : 1;
    bh_duty *= f_mass;
    if (bh_duty < 1e-4) bh_duty = 1e-4;

    //fprintf(stderr, "mh=%f, sm=%f, obs_sm=%f, mbh=%f\n", logm, sm, obs_sm, mbh);  

    
    float mf = (logm - M_MIN)*BPDEX;
    int64_t j = mf;
    mf -= j;
    float sfr = biterp(steps[ts].sfr[j], steps[ts+1].sfr[j],
       steps[ts].sfr[j+1], steps[ts+1].sfr[j+1],
       tf, mf);
    float scatter_corr = log10(exp(pow(ins_scatter*log(10), 2)/2.0));
    sfr -= scatter_corr;
    sfr += ins_scatter*(sm_offset*rho + sfr_offset*sqrt(1.0-rho*rho));
    if (halos[i].upid >= 0) { sfr += log10(halos[i].m/halos[i].mp); }
    // if (sfr - sm > 11) obs_sm += kappa;
    float obs_sfr = sfr + mu + kappa;

    // Calculate the cumulative probabilities of hosting AGN above certain Eddington ratio
    // or luminosity thresholds.
    for (int k=0; k<5; k++)
    {
      // Eddington ratio thresholds. Only calculate when k < 4, because we only have
      // 4 threshold values.
      
      float eta_crit;
      double prob1, prob2, bher_f;
      int64_t bher_b;

      float mf = (logm - M_MIN) * BPDEX;
      int64_t mb = mf;
      mf -= mb;
      if (mb == M_BINS - 1) {mb = M_BINS - 2; mf = 1;}


      if (k < 4)
      {

        eta_crit = -3.0 + k;

        // for the (mb)-th halo mass bin
        prob1 = 0.0;
        bher_f = (eta_crit - steps[ts].bh_eta[mb] - steps[ts].ledd_min[mb]) * steps[ts].ledd_bpdex[mb];
        bher_b = bher_f;
        bher_f -= bher_b;
        if (bher_b >= BHER_BINS - 1) {bher_b = BHER_BINS - 2; bher_f = 1;}

        prob1 += (1 - bher_f) * steps[ts].bher_dist[mb*BHER_BINS+bher_b];
        for (int ii=bher_b+1; ii<BHER_BINS; ii++)
          prob1 += steps[ts].bher_dist[mb*BHER_BINS+ii];
        if (steps[ts].bher_dist_norm[mb] > 0) prob1 /= steps[ts].bher_dist_norm[mb];

        // for the (mb+1)-th halo mass bin
        prob2 = 0.0;
        bher_f = (eta_crit - steps[ts].bh_eta[mb+1] - steps[ts].ledd_min[mb+1]) * steps[ts].ledd_bpdex[mb+1];
        bher_b = bher_f;
        bher_f -= bher_b;
        if (bher_b >= BHER_BINS - 1) {bher_b = BHER_BINS - 2; bher_f = 1;}

        prob2 += (1 - bher_f) * steps[ts].bher_dist[(mb+1)*BHER_BINS+bher_b];
        for (int ii=bher_b+1; ii<BHER_BINS; ii++)
          prob2 += steps[ts].bher_dist[(mb+1)*BHER_BINS+ii];
        if (steps[ts].bher_dist_norm[(mb+1)] > 0) prob2 /= steps[ts].bher_dist_norm[(mb+1)];

        // interpolate between the mb-th and the (mb+1)-th mass bin.
        p_eta[k] = (prob1 + mf * (prob2 - prob1)) * bh_duty;
      }
      

      //fprintf(stderr, "mbh=%f, eta_crit=%f, p_eta[%d]=%e, prob1=%e, prob2=%e, bher_b=%d, bher_f=%f, bher_dist_norm=%e\n", mbh, eta_crit, k, p_eta[k], prob1, prob2, bher_b, bher_f, steps[ts].bher_dist_norm[mb+1]); 

      // luminosity thresholds.
      float lbol_crit = 41.0 + k;
      eta_crit = lbol_crit - 38.1 - mbh;


      // for the mb-th halo mass bin
      prob1 = 0.0;
      bher_f = (eta_crit - steps[ts].bh_eta[mb] - steps[ts].ledd_min[mb]) * steps[ts].ledd_bpdex[mb];
      bher_b = bher_f;
      bher_f -= bher_b;
      if (bher_b >= BHER_BINS - 1) {bher_b = BHER_BINS - 2; bher_f = 1;}

      prob1 += (1 - bher_f) * steps[ts].bher_dist[mb*BHER_BINS+bher_b];
      for (int ii=bher_b+1; ii<BHER_BINS; ii++)
        prob1 += steps[ts].bher_dist[mb*BHER_BINS+ii];
      if (steps[ts].bher_dist_norm[mb] > 0) prob1 /= steps[ts].bher_dist_norm[mb];

      // for the (mb+1)-th halo mass bin
      prob2 = 0.0;
      bher_f = (eta_crit - steps[ts].bh_eta[mb+1] - steps[ts].ledd_min[mb+1]) * steps[ts].ledd_bpdex[mb+1];
      bher_b = bher_f;
      bher_f -= bher_b;
      if (bher_b >= BHER_BINS - 1) {bher_b = BHER_BINS - 2; bher_f = 1;}

      prob2 += (1 - bher_f) * steps[ts].bher_dist[(mb+1)*BHER_BINS+bher_b];
      for (int ii=bher_b+1; ii<BHER_BINS; ii++)
        prob2 += steps[ts].bher_dist[(mb+1)*BHER_BINS+ii];
      if (steps[ts].bher_dist_norm[(mb+1)] > 0) prob2 /= steps[ts].bher_dist_norm[(mb+1)];

      // interpolate between the mb-th and the (mb+1)-th mass bin.
      p_lbol[k] = (prob1 + mf * (prob2 - prob1)) * bh_duty;
    }

    float eta_lbol_Aird = lbol_aird - 38.1 - mbh;
    // If the Eddington ratio that gives logLbol = 45 is bigger than 0.1, then we 
    // take the prob value for logLbol = 45.
    p_aird = eta_lbol_Aird > eta_aird ? p_lbol[4] : p_eta[2];

    fprintf(output, "%.6f %8"PRId64" %8"PRId64" %8"PRId64" %d %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6e %6e %.6e %.6e %.6e %.6f %.6f %.6f %.6e %.6e %.6e %.6e %.6e %.6f %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
	    scale, halos[i].id, halos[i].descid, halos[i].upid, halos[i].flags,
      halos[i].uparent_dist,
      halos[i].pos[0], halos[i].pos[1], halos[i].pos[2], 
      halos[i].pos[3], halos[i].pos[4], halos[i].pos[5], 
      halos[i].m, halos[i].v, halos[i].mp, halos[i].vmp,
      halos[i].r, halos[i].rank1, halos[i].rank2, A_UV(halos[i]),
      exp10(sm), halos[i].icl, exp10(sfr),
      exp10(obs_sm), exp10(obs_sfr), halos[i].obs_uv, exp10(mbh),
      p_eta[0], p_eta[1], p_eta[2], p_eta[3], p_lbol[0], p_lbol[1], p_lbol[2], p_lbol[3], p_lbol[4], p_aird);
  }
  fclose(output);
}


int main(int argc, char **argv)
{
  int i;
  enum parsetype t[NUM_PARAMS];
  void *data[NUM_PARAMS];
  char buffer[2048] = {0};
  struct smf_fit base_smf;

  if (argc<3) {
    fprintf(stderr, "Usage: %s mass_cache catalog1 ... < mcmc_output\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++) {
    t[i] = PARSE_FLOAT64;
    data[i] = &(base_smf.params[i]);
    // printf("data[%d]=%.6f\n", i, data[i]);
  }

  fgets(buffer, 2048, stdin);
  if (stringparse(buffer, data, t, NUM_PARAMS) != NUM_PARAMS) {
    fprintf(stderr, "%d\n", stringparse(buffer, data, t, NUM_PARAMS));
    fprintf(stderr, "Could not parse input!\n");
    exit(1);
  }

  r250_init(87L);
  gen_exp10cache();
  setup_psf(1);
  gsl_set_error_handler_off();
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(base_smf) = 0;
  calc_sfh(&base_smf);

#pragma omp for schedule(guided,5)
  for (i=2; i<argc; i++)
    add_agn_to_catalog(base_smf, argv[i]);
  return 0;
}
