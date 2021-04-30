#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>
#include "check_syscalls.h"
#include "make_sf_catalog.h"

#define H0 0.7

int main(int argc, char **argv) {
  int64_t i, new_length;
  if (argc<2) {
    printf("Usage: %s acc_list.bin\n", argv[0]);
    exit(1);
  }

  struct catalog_halo *p = check_mmap_file(argv[1], 'r', &new_length);
  int64_t num_p = new_length / sizeof(struct catalog_halo);
  printf("# Column 1: ID\n");
  printf("# Column 2: sat?\n");
  printf("# Column 3-5: x,y,z\n");
  printf("# Column 4-6: Vx,Vy,Vz\n");
  printf("# Column 7: Vpeak\n");
  printf("# Column 8: Ignore\n");
  printf("# Column 9: log10(M*)\n");
  printf("# Column 10: log10(sSFR)\n");
  printf("# Column 11: iquench (=1 for quenched)\n");

  double sm_thresh = pow(10, 9.8);

  for (i=0; i<num_p; i++) {
    struct catalog_halo *tp = p+i;
    if (tp->sm < sm_thresh) continue;
    printf("%"PRId64" %f %f %f %f %f %f %f %f 0 %f %f %d\n",
	   tp->id, tp->sat, tp->pos[0], tp->pos[1], tp->pos[2], tp->pos[3], tp->pos[4], tp->pos[5], tp->vp, log10(tp->sm), log10(tp->obs_sfr/tp->sm), (tp->obs_sfr < 1e-11*tp->sm) ? 1 : 0);
  }
  munmap(p, new_length);
  return 0;
}
