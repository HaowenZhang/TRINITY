#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "check_syscalls.h"
#include "mt_rand.h"

#define NUM_DRAWS 22

struct pair {
  float sm, re, prob;
};

struct pair *pairs = NULL;
int64_t num_pairs = 0;

struct pair draws[NUM_DRAWS];

void draw_set(void) {
  int64_t i;
  while (i<NUM_DRAWS) {
    int64_t index = num_pairs*dr250();
    float prob = dr250();
    if (prob > pairs[index].prob) continue;
    draws[i] = pairs[index];
    i++;
  }
}

float kendalls_tau(void) {
  int64_t i,j;
  draw_set();
  int64_t cp=0, dp=0;
  double c = 0;
  for (i=0; i<NUM_DRAWS; i++) {
    for (j=i+1; j<NUM_DRAWS; j++) {
      c++;
      if ((draws[i].sm-draws[j].sm)*(draws[i].re-draws[j].re)>0) cp++;
      else dp++;
    }
  }
  return ((cp-dp)/c);
}


int main(void) {
  struct pair p = {0};
  FILE *input = check_fopen("grb_sm_offsets_filtered.dat", "r");
  char buffer[1024];
  r250_init(87L);
  float max_prob = 0;

  while (fgets(buffer, 1024, input)) {
    if (sscanf(buffer, "%f %f %f", &p.sm, &p.re, &p.prob)<3) continue;
    if (!(num_pairs%1000))
      pairs = check_realloc(pairs, sizeof(struct pair)*(num_pairs+1000), "");
    pairs[num_pairs] = p;
    num_pairs++;
    if (p.prob > max_prob) max_prob = p.prob;
  }
  fclose(input);

  int64_t i;
  for (i=0; i<num_pairs; i++) pairs[i].prob /= max_prob;
  for (i=0; i<10000; i++) printf("%f\n", kendalls_tau());
  return 0;
}
