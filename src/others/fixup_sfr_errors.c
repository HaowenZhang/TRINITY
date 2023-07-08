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

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

float m_to_bin(float m) {
  m -= INV_BPDEX/2.0;
  if (m<M_MIN) return 0;
  if (m>(M_MAX-INV_BPDEX)) return (M_BINS-2);
  return ((m-M_MIN)*BPDEX);
}

float sm_at_m(float m, int i) {
  float f = m_to_bin(m);
  int b = f;
  f-=b;
  if (!b) return 0;
  float sm = steps[i].sm[b] + f*(steps[i].sm[b+1]-steps[i].sm[b]);
  return sm;
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

#define NUM_MS 5

static void read_params(char *buffer, double *data, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    data[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}

int main(int argc, char **argv)
{
  float m;
  double fb = 0.0455/(0.0455+0.226);
  char buffer[1024];
  struct smf_fit the_smf;
  int64_t type;
  FILE *input;

  if (argc<6) {
    fprintf(stderr, "Usage: %s mass mass_cache type stats_output file\n", argv[0]);
    exit(1);
  }

  m = atof(argv[1]);
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();

  input = check_fopen(argv[4], "r");
  fgets(buffer, 1024, input);
  fclose(input);
  read_params(buffer, the_smf.params, NUM_PARAMS+2);
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  type = atol(argv[3]);
  input = check_fopen(argv[5], "r");
  while (fgets(buffer, 1024, input)) {
    double scale, sfr, up, down, f;
    sscanf(buffer, "%lf %lf %lf %lf\n", &scale, &sfr, &up, &down);
    int64_t step;
    calc_step_at_z(1.0/scale-1.0, &step, &f);
    step += f+0.5;
    double sfr_err=0;
    double rel_sfr_err = 2.5e-11;

    if (type==0) { //Plain Sfr
      double sm = calc_sm_at_m(m, steps[step].smhm);
      sfr_err = pow(10,sm)*rel_sfr_err;
    }

    else if (type==1) { //SFH
      double sm = calc_sm_at_m(m_at_a(m, steps[step].scale), steps[step].smhm);
      sfr_err = pow(10,sm)*rel_sfr_err;
    }
    
    else if (type==2) { //SFR / MAR / fb (Mnow)
      double sm = calc_sm_at_m(m, steps[step].smhm);
      sfr_err = pow(10,sm)*rel_sfr_err / mar_from_mbins(step, (m-M_MIN)*BPDEX) / fb;
    }

    else if (type==3) { //SFR / MAR / fb (Mthen)
      double sm = calc_sm_at_m(m_at_a(m, steps[step].scale), steps[step].smhm);
      sfr_err = pow(10,sm)*rel_sfr_err / ma_rate(m, steps[step].scale) / fb;
    }

    sfr_err = 0;

    if (down < sfr_err) down = sfr_err;
    if (up < sfr_err) up = sfr_err;
    if (down > 0.999*sfr) down = 0.999*sfr;
    printf("%f %f %f %f\n", scale, sfr, up, down);
  }
  fclose(input);
  return 0;
}
