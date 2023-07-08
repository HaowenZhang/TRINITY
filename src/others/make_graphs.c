#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "check_syscalls.h"
#include "make_graphs.h"
#include "mah.h"

extern int64_t num_outputs;
extern struct timestep *steps;

int main(int argc, char **argv)
{
  int64_t i;
  char buffer[1024];
  FILE *ifile;
  struct smf_fit input;

  if (argc<3) {
    fprintf(stderr, "Usage: %s mass_cache mcmc_file\n", argv[0]);
    exit(1);
  }

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  ifile = check_fopen(argv[2], "r");

  i=0;
  while (fgets(buffer, 1024, ifile)) {
    input = parse_smf_fit(buffer);
    calc_sfh(&input);
    print_sfh(i++, input);
    if (i==1000) break;
  }
  fclose(ifile);
  return 0;
}


#define B_START 6
#define B_END 12
#define B_BPDEX 10
#define B_NB ((B_END-B_START)*B_BPDEX+1)
float conv_coeffs[B_NB][M_BINS];

void gen_conv_coeffs(struct smf_fit f) {
  int64_t i,j;
  struct smf smf = smhm_at_z(1.0/steps[num_outputs-1].scale-1.0, f);
  double scatter2 = smf.scatter*smf.scatter + smf.obs_scatter*smf.obs_scatter;
  for (i=0; i<M_BINS; i++) for (j=0; j<B_NB; j++) conv_coeffs[j][i] = 0;
  for (i=0; i<M_BINS; i++) {
    if (steps[num_outputs-1].sm[i] <= 1) continue;
    if (i*INV_BPDEX + M_MIN > 15) continue;
    double log_hsm = log10(steps[num_outputs-1].sm[i]);
    for (j=0; j<B_NB; j++) {
      double sm = B_START + j/((double)B_BPDEX);
      double delta_m = sm - log_hsm;
      double psf = exp(-0.5*(delta_m*delta_m/scatter2));
      //New Kappa:
      double passive_frac = 1.0/(exp10(-1.3*(sm-smf.passive_mass))+1.0);
      delta_m -= smf.kappa;
      psf = psf*passive_frac + 
	(1.0-passive_frac)*exp(-0.5*(delta_m*delta_m/scatter2));
      conv_coeffs[j][i] = psf;
    }
  }
  for (j=0; j<B_NB; j++) {
    double sum = 0;
    for (i=0; i<M_BINS; i++) sum += conv_coeffs[j][i];
    if (sum>0)
      for (i=0; i<M_BINS; i++) conv_coeffs[j][i] /= sum;
  }
}

float find_fraction_dt(int i, float f, float total_sm) {
  int64_t j=0, j2=0;
  double sum = 0, min_dt;
  for (j2=0; j2<num_outputs; j2++) {
    sum += steps[num_outputs-1].sm_hist[i*num_outputs+j2]*steps[num_outputs-1].smloss[j2];
    if (sum > f*(total_sm)) break;
  }
    min_dt = (j2-j)*steps[num_outputs-1].dt*(f*(total_sm)/sum);
  for (j=0; j<num_outputs; j++) {
    sum -= steps[num_outputs-1].sm_hist[i*num_outputs+j]*steps[num_outputs-1].smloss[j];
    sum -= steps[num_outputs-1].sm_hist[i*num_outputs+j2]*steps[num_outputs-1].smloss[j2];
    for (; j2<num_outputs; j2++) {
      sum += steps[num_outputs-1].sm_hist[i*num_outputs+j2]*steps[num_outputs-1].smloss[j2];
      if (sum > f*(total_sm)) break;
    }
    if (j2==num_outputs) break;
    double dt = (j2-j)*steps[num_outputs-1].dt*(f*(total_sm)/sum);
    if (dt < min_dt) {
      min_dt = dt;
      /*      if (min_dt < 1.2e9) {
	printf("Mass: %f; Scale 1: %f; Scale 2: %f; Sum: %e; Total: %e\n", M_MIN + i*INV_BPDEX, steps[j].scale, steps[j2].scale, sum, total_sm);
	}*/
    }
  }
  return min_dt;
}

double convolve(float *array, int64_t j) {
  int64_t i;
  double sum=0;
  for (i=0; i<M_BINS; i++) sum += array[i]*conv_coeffs[j][i];
  return sum;
}

float m_to_bin(float m) {
  m -= INV_BPDEX/2.0;
  if (m<M_MIN) return 0;
  if (m>(M_MAX-INV_BPDEX)) return (M_BINS-2);
  return ((m-M_MIN)*BPDEX);
}

float _mar_from_mbins(int64_t n, int64_t j) {
  int64_t i;
  if (!n) return pow(10, M_MIN+(j+0.5)*INV_BPDEX)/steps[n].dt;
  if (n>=num_outputs-1) return _mar_from_mbins(num_outputs-2, j);
  if (j>=M_BINS-1) return 0;
  //Find out average progenitor mass:
  double sum = 0;
  double count = 0;
  for (i=0; i<M_BINS; i++) {
    if (!steps[n].mmp[i*M_BINS + j]) continue;
    sum += pow(10, M_MIN+(i+0.5)*INV_BPDEX)*steps[n].mmp[i*M_BINS+j];
    count += steps[n].mmp[i*M_BINS+j];
  }
  if (!count) return 0;
  sum /= count;
  return ((pow(10, M_MIN+(j+0.5)*INV_BPDEX) - sum)/steps[n].dt);
}

float mar_from_mbins(int64_t n, int64_t j) {
  float mar1 = _mar_from_mbins(n,j);
  float mar2 = _mar_from_mbins(n+1,j);
  return (0.5*(mar1+mar2));
}



void print_sfh(int64_t count, struct smf_fit input) {
  char buffer[1024];
  int64_t i,j,k,step;
  static FILE *output[50];
  double output_zs[] = {8,7,6,5,4,3,2,1,0.5,0.001};
  int num_z_outputs = 10;
  double scatter_corr = exp(pow(SCATTER(input)*log(10), 2)/2.0);
  double total_sm[M_BINS];

  sprintf(buffer, "results/smf_fit");
  if (!count) {
    output[0] = check_fopen(buffer, "w");
    for (i=0; i<NUM_PARAMS; i++) fprintf(output[0], "%.12f ", input.params[i]);
    fclose(output[0]);
  }

  sprintf(buffer, "results/sfr.dat");
  if (!count) output[0] = check_fopen(buffer, "w");
  fprintf(output[0], "#Scale ");
  for (j=0; j<M_BINS; j++)
    fprintf(output[0], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[0], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[0], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j++) {
      fprintf(output[0], "%e ", steps[i].sfr[j]);
    }
    fprintf(output[0], "\n");
  }

  sprintf(buffer, "results/sm.dat");
  if (!count) output[1] = check_fopen(buffer, "w");
  fprintf(output[1], "#Scale ");
  for (j=0; j<M_BINS; j+=BPDEX)
    fprintf(output[1], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[1], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[1], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j+=BPDEX) {
      fprintf(output[1], "%e ", steps[i].sm[j]);
    }
    fprintf(output[1], "\n");
  }

  sprintf(buffer, "results/sm_full.dat");
  if (!count) output[20] = check_fopen(buffer, "w");
  fprintf(output[20], "#Scale ");
  for (j=0; j<M_BINS; j++)
    fprintf(output[20], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[20], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[20], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j++) {
      fprintf(output[20], "%e ", steps[i].sm[j]);
    }
    fprintf(output[20], "\n");
  }

  sprintf(buffer, "results/smhist.dat");
  if (!count) output[2] = check_fopen(buffer, "w");
  fprintf(output[2], "#Scale ");
  for (j=0; j<M_BINS; j++)
    fprintf(output[2], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[2], "\n");
  memset(total_sm, 0, sizeof(double)*M_BINS);
  for (i=0; i<num_outputs; i++) {
    fprintf(output[2], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j++) {
      total_sm[j] += steps[num_outputs-1].smloss[i]*steps[num_outputs-1].sm_hist[j*num_outputs+i];
      fprintf(output[2], "%e ", total_sm[j]);
    }
    fprintf(output[2], "\n");
  }

  sprintf(buffer, "results/icl.dat");
  if (!count) output[3] = check_fopen(buffer, "w");
  fprintf(output[3], "#Scale ");
  for (j=0; j<M_BINS; j++)
    fprintf(output[3], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[3], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[3], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j++) {
      fprintf(output[3], "%e ", steps[i].sm_icl[j]);
    }
    fprintf(output[3], "\n");
  }

  sprintf(buffer, "results/icl_prog.dat");
  if (!count) output[4] = check_fopen(buffer, "w");
  fprintf(output[4], "#Scale ");
  for (j=0; j<M_BINS; j++)
    fprintf(output[4], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[4], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[4], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j++) {
      float m = m_at_a(M_MIN+j*INV_BPDEX, steps[i].scale);
      if (m>(M_MAX-2*INV_BPDEX)) m=M_MAX-(2*INV_BPDEX);
      float f = (m-M_MIN)*BPDEX;
      if (f<0) f=0;
      int64_t j2 = f;
      f-=j2;
      fprintf(output[4], "%e ", steps[i].sm_icl[j2] + f*(steps[i].sm_icl[j2+1]-steps[i].sm_icl[j2]));
    }
    fprintf(output[4], "\n");
  }


  sprintf(buffer, "results/smhm.dat");
  if (!count) output[5] = check_fopen(buffer, "w");
  fprintf(output[5], "#Mass Smhm\n");
  i=num_outputs-1;
  for (j=0; j<M_BINS; j++) {
    double hm = M_MIN+(j+0.5)*INV_BPDEX;
    fprintf(output[5], "%e %e\n", hm, log10(steps[i].sm[j])-hm);
  }

  sprintf(buffer, "results/madau.dat");
  if (!count) output[6] = check_fopen(buffer, "w");
  fprintf(output[6], "#Redshift CSFR\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[6], "%f ", 1.0/steps[i].scale-1.0);
    fprintf(output[6], "%e\n", steps[i].observed_cosmic_sfr);
  }

  sprintf(buffer, "results/cosmic_sfr.dat");
  if (!count) output[7] = check_fopen(buffer, "w");
  fprintf(output[7], "#Redshift");
  for (j=0; j<M_BINS; j++)
    fprintf(output[7], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[7], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[7], "%f ", 1.0/steps[i].scale-1.0);
    for (j=0; j<M_BINS; j++) {
      fprintf(output[7], "%e ", steps[i].sfr[j]*steps[i].t[j]/scatter_corr);
    }
    fprintf(output[7], "\n");
  }

  sprintf(buffer, "results/cosmic_indv_sm.dat");
  if (!count) output[8] = check_fopen(buffer, "w");
  fprintf(output[8], "#Redshift");
  for (j=0; j<M_BINS; j++)
    fprintf(output[8], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[8], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[8], "%f ", 1.0/steps[i].scale-1.0);
    for (j=0; j<M_BINS; j++) {
      fprintf(output[8], "%e ", steps[i].sm[j]*steps[i].t[j]);
    }
    fprintf(output[8], "\n");
  }

  sprintf(buffer, "results/cosmic_sm.dat");
  if (!count) output[9] = check_fopen(buffer, "w");
  fprintf(output[9], "#Redshift Mass\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[9], "%f ", 1.0/steps[i].scale-1.0);
    double smt = 0;
    for (j=0; j<M_BINS; j++)
      smt+=steps[i].sm[j]*steps[i].t[j];
    fprintf(output[9], "%e\n", smt*scatter_corr);
  }

  sprintf(buffer, "results/ssfr.dat");
  if (!count) output[10] = check_fopen(buffer, "w");
  fprintf(output[10], "#Scale ");
  for (j=0; j<M_BINS; j+=BPDEX)
    fprintf(output[10], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[10], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[10], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j+=BPDEX) {
      double ssfr = (steps[i].sm[j] > 0) ? (steps[i].sfr[j] / steps[i].sm[j]) : 0;
      fprintf(output[10], "%e ", ssfr);
    }
    fprintf(output[10], "\n");
  }

  sprintf(buffer, "results/ssfr_sm.dat");
  {
   double sm_lim = 12;
   double sm_start = 8;
   double sm, ssfr;
   if (!count) output[11] = check_fopen(buffer, "w");
   fprintf(output[11], "#Scale ");
   for (j = (sm_start - SM_MIN)*SM_BPDEX; ; j+=10) {
     sm = SM_MIN + j*SM_INV_BPDEX;
     if ((sm > sm_lim) || (j >= SM_BINS + SM_EXTRA)) break;
     fprintf(output[11], "%.1f ", sm);
   }
   fprintf(output[11], "\n");

   for (i=0; i<num_outputs; i++) {
     fprintf(output[11], "%f ", steps[i].scale);
     for (j = (sm_start - SM_MIN)*SM_BPDEX; ; j+=10) {
       sm = SM_MIN + j*SM_INV_BPDEX;
       if ((sm > sm_lim) || (j >= SM_BINS + SM_EXTRA)) break;
       ssfr = steps[i].sfr_sm[j+SM_EXTRA]/pow(10,sm);
       fprintf(output[11], "%e ", ssfr);
     }
     fprintf(output[11], "\n");
   }
 }

  sprintf(buffer, "results/ages.dat");
  if (!count) output[12] = check_fopen(buffer, "w");
  float age;
  float ages[M_BINS], time50[M_BINS], time90[M_BINS];
  gen_conv_coeffs(input);
  for (i=0; i<M_BINS; i++) {
    float dt = 0;
    float total_sm = 0;
    if (i*INV_BPDEX + M_MIN > 15) continue;
    age = 0;
    if (steps[num_outputs-1].sm[i] <= 0) continue;
    for (j=num_outputs-1; j>=0; j--) {
      dt += steps[j].dt;
      age += steps[num_outputs-1].sm_hist[i*num_outputs+j]*steps[num_outputs-1].smloss[j]
	*(dt);
      total_sm += steps[num_outputs-1].sm_hist[i*num_outputs+j]*steps[num_outputs-1].smloss[j];
    }
    age /= total_sm;
    ages[i] = age;
    time50[i] = find_fraction_dt(i, 0.5, total_sm);
    time90[i] = find_fraction_dt(i, 0.9, total_sm);
  }
  for (i=0; i<B_NB; i++) {
    double sm = pow(10, B_START + i/((double)B_BPDEX));
    fprintf(output[12], "%e %g %g %g\n", sm, convolve(ages, i), convolve(time90, i), convolve(time50,i));
  }

  sprintf(buffer, "results/sfr_ma_atz.dat");
  if (!count) output[13] = check_fopen(buffer, "w");
  fprintf(output[13], "#Scale ");
  for (j=0; j<M_BINS; j+=BPDEX)
    fprintf(output[13], "%.1f ", M_MIN+j*INV_BPDEX);
  fprintf(output[13], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[13], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j+=BPDEX) {
      double fb = 0.0455/(0.0455+0.226);
      double sfr_ma_atz = steps[i].sfr[j]/mar_from_mbins(i,j)/fb;
      fprintf(output[13], "%e ", sfr_ma_atz);
    }
    fprintf(output[13], "\n");
  }  

  sprintf(buffer, "results/sfr_ma_atz0.dat");
  if (!count) output[14] = check_fopen(buffer, "w");
  sprintf(buffer, "results/smhm_hist.dat");
  if (!count) output[15] = check_fopen(buffer, "w");
  sprintf(buffer, "results/sm_hist_rel.dat");
  if (!count) output[16] = check_fopen(buffer, "w");
  fprintf(output[14], "#Scale ");
  fprintf(output[15], "#Scale ");
  fprintf(output[16], "#Scale ");
  for (j=0; j<M_BINS; j+=BPDEX) {
    fprintf(output[14], "%.1f ", M_MIN+j*INV_BPDEX);
    fprintf(output[15], "%.1f ", M_MIN+j*INV_BPDEX);
    fprintf(output[16], "%.1f ", M_MIN+j*INV_BPDEX);
  }
  fprintf(output[14], "\n");
  fprintf(output[15], "\n");
  fprintf(output[16], "\n");
  for (i=0; i<num_outputs; i++) {
    fprintf(output[14], "%f ", steps[i].scale);
    fprintf(output[15], "%f ", steps[i].scale);
    fprintf(output[16], "%f ", steps[i].scale);
    for (j=0; j<M_BINS; j+=BPDEX) {
      double m = M_MIN+j*INV_BPDEX;
      double m_then = m_at_a(m, steps[i].scale);
      double f = m_to_bin(m_then);
      float mar = ma_rate(m, steps[i].scale);
      float sm = pow(10, calc_sm_at_m(m_then, steps[i].smhm));
      int64_t b = f;
      f-=b;
      double sfr = steps[i].sfr[b] + f*(steps[i].sfr[b+1]-steps[i].sfr[b]);
      double fb = 0.0455/(0.0455+0.226);
      double sfr_ma_atz0 = sfr/mar/fb;
      fprintf(output[14], "%e ", sfr_ma_atz0);
      fprintf(output[15], "%e ", sm/pow(10, m_then));
      fprintf(output[16], "%e ", sm/pow(10,calc_sm_at_m(m, steps[num_outputs-1].smhm)));
    }
    fprintf(output[14], "\n");
    fprintf(output[15], "\n");
    fprintf(output[16], "\n");
  }  
  

  sprintf(buffer, "results/ssfr_sm_obs.dat");
  {
   double sm_lim = 12;
   double sm_start = 8;
   double sm, ssfr;
   if (!count) output[17] = check_fopen(buffer, "w");
   fprintf(output[17], "#Scale ");
   for (j = (sm_start - SM_MIN)*SM_BPDEX; ; j+=10) {
     sm = SM_MIN + j*SM_INV_BPDEX;
     if ((sm > sm_lim) || (j >= SM_BINS + SM_EXTRA)) break;
     fprintf(output[17], "%.1f ", sm);
   }
   fprintf(output[17], "\n");

   for (i=0; i<num_outputs; i++) {
     fprintf(output[17], "%f ", steps[i].scale);
     for (j = (sm_start - SM_MIN)*SM_BPDEX; ; j+=10) {
       sm = SM_MIN + j*SM_INV_BPDEX;
       if ((sm > sm_lim) || (j >= SM_BINS + SM_EXTRA)) break;
       ssfr = calc_ssfr(sm, 1.0/steps[i].scale-1.0);
       fprintf(output[17], "%e ", ssfr);
     }
     fprintf(output[17], "\n");
   }
 }

  sprintf(buffer, "results/sfr_ma_time_weighted.dat");
  if (!count) output[18] = check_fopen(buffer, "w");
  fprintf(output[18], "#Mass SFR_MA\n");
  for (j=0; j<M_BINS; j++) {
    float avg_sfr_ma = 0, total = 0;
    for (i=0; i<num_outputs; i++) {
      if (steps[i].scale < 0.2) continue;
      double fb = 0.0455/(0.0455+0.226);
      double sfr_ma_atz = steps[i].sfr[j]/mar_from_mbins(i,j)/fb;
      if (sfr_ma_atz > 0 && sfr_ma_atz < 1) {
	total++;
	avg_sfr_ma += sfr_ma_atz;
      }
    }
    if (total)
      fprintf(output[18], "%f %e\n", (j+0.5)*INV_BPDEX+M_MIN, avg_sfr_ma/total);
  }


 sprintf(buffer, "results/sfr_ma_csfr_weighted.dat");
  if (!count) output[19] = check_fopen(buffer, "w");
  fprintf(output[19], "#Mass SFR_MA\n");
  for (j=0; j<M_BINS; j++) {
    float avg_sfr_ma = 0, total = 0;
    for (i=0; i<num_outputs; i++) {
      if (steps[i].scale < 0.2) continue;
      double fb = 0.0455/(0.0455+0.226);
      double sfr_ma_atz = steps[i].sfr[j]/mar_from_mbins(i,j)/fb;
      if (sfr_ma_atz > 0 && sfr_ma_atz < 1) {
	total+=steps[i].cosmic_sfr;
	avg_sfr_ma += sfr_ma_atz*steps[i].cosmic_sfr;
      }
    }
    if (total)
      fprintf(output[19], "%f %e\n", (j+0.5)*INV_BPDEX+M_MIN, avg_sfr_ma/total);
  }


  int64_t fd_start = 21;
  for (k=0; k<num_z_outputs; k++) {
    sprintf(buffer, "results/sfh_z%.1f.dat", output_zs[k]);
    if (!count) output[fd_start+k] = check_fopen(buffer, "w");
    step = find_step(1.0/(output_zs[k]+1.0));

    fprintf(output[fd_start+k], "#Scale ");
    for (j=0; j<M_BINS; j++)
      fprintf(output[fd_start+k], "%.1f ", M_MIN+j*INV_BPDEX);
    fprintf(output[fd_start+k], "\n");

    for (i=0; i<=step; i++) {
      fprintf(output[fd_start+k], "%f ", steps[i].scale);
      for (j=0; j<M_BINS; j++) {
	fprintf(output[fd_start+k], "%e ", steps[step].sm_hist[j*num_outputs+i] / steps[i].dt);
      }
      fprintf(output[fd_start+k], "\n");
    }
  }



  if (count==999) {
    for (i=0; i<fd_start+num_z_outputs; i++) fclose(output[i]);
  }
}

int64_t find_step(double scale) {
  int64_t i;
  for (i=0; i<num_outputs; i++)
    if (steps[i].scale>=scale) break;

  return i;
}


static void read_params(char *buffer, double *data, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  double val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    data[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}

struct smf_fit parse_smf_fit(char *buffer) {
  struct smf_fit a;
  read_params(buffer, a.params, NUM_PARAMS);
  return a;
}


