#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include "make_sf_catalog.h"
#include "check_syscalls.h"
#include "stringparse.h"

struct particle {
  float pos[2];
};

#define FAST3TREE_DIM 2
#define FAST3TREE_TYPE struct particle
#include "fast3tree.c"

struct particle *p = NULL;
int64_t num_p = 0;

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s box_size particles acc_list\n", argv[0]);
    exit(1);
  }

  float box_size = atof(argv[1]);
  char buffer[1024];
  FILE *in = check_fopen(argv[2], "r");
  struct particle tp;
  void *data[2] = {tp.pos, tp.pos+1};
  enum parsetype pt[2] = {PARSE_FLOAT32, PARSE_FLOAT32};
  while (fgets(buffer, 1024, in)) {
    if (stringparse(buffer, data, pt, 2)!=2) continue;
    check_realloc_every(p, sizeof(struct particle), num_p, 1000);
    p[num_p] = tp;
    num_p++;
  }
  fclose(in);

  struct fast3tree *tree = fast3tree_init(num_p, p);
  struct fast3tree_results *res = fast3tree_results_init();
  _fast3tree_set_minmax(tree, 0, box_size);

  int64_t i, length = 0;
  struct packed_halo *ph = check_mmap_file(argv[3], 'r', &length);
  mlock(ph, length);
  assert(!(length % sizeof(struct packed_halo)));
  int64_t num_halos = length / sizeof(struct packed_halo);

#define BPDEX 8
#define MIN_R -2
#define NUM_BINS 20
  printf("#ID X Y Z Particle Counts\n#Bin starts: 0");
  for (i=0; i<NUM_BINS-1; i++)
    printf(" %f", pow(10, MIN_R+i/(double)BPDEX));
  printf("\n#Bin ends:");
  for (i=0; i<NUM_BINS; i++)
    printf(" %f", pow(10, MIN_R+i/(double)BPDEX));
  printf("\n");
  printf("#All units: comoving Mpc/h\n");

  for (i=0; i<num_halos; i++) {
    fast3tree_find_sphere_periodic(tree, res, ph[i].pos, pow(10, MIN_R+(NUM_BINS-1)/(double)BPDEX));
    int64_t j;
    int64_t bins[NUM_BINS] = {0};
    printf("%"PRId64" %f %f %f", ph[i].id, ph[i].pos[0], ph[i].pos[1], ph[i].pos[2]);
    int64_t k=0;
    for (j=0; j<res->num_points; j++) {
      double dx, ds=0;
      for (k=0; k<2; k++) {
	dx = fabs(ph[i].pos[k] - res->points[j]->pos[k]);
	if (dx > box_size/2.0) dx = box_size - dx;
	ds += dx*dx;
      }
      int64_t bin = (0.5*log10(ds)-MIN_R)*(BPDEX)+1;
      if (bin<0) bin = 0;
      assert(bin < NUM_BINS);
      bins[bin]++;
    }
    for (k=0; k<NUM_BINS; k++)
      printf(" %"PRId64, bins[k]);
    printf("\n");
  }
  munmap(ph, length);
  return 0;
}
