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
    printf("Usage: %s am_catalog.dat\n", argv[0]);
    exit(1);
  }

  FILE *in = check_fopen(argv[1], "r");
  FILE *out = check_fopen("am_catalog.bin", "w");
  struct catalog_halo p = {0};
  char buffer[1024];
  while (fgets(buffer, 1024, in)) {
    int64_t upid;
    if (buffer[0] == '#') continue;
    float sm, ssfr;
    if (sscanf(buffer, "%"SCNd64" %"SCNd64" %f %f %f %f %f %f %f %f %f %f",
	       &p.id, &upid, p.pos, p.pos+1, p.pos+2, p.pos+3, p.pos+4, p.pos+5, 
	       &p.vp, &p.mp, &sm, &ssfr) != 12) continue;
    p.sm = pow(10, sm);
    p.obs_sfr = p.sfr = pow(10, sm+ssfr);
    p.sat = (upid < 0) ? 0.0 : 1.0;
    fwrite(&p, sizeof(struct catalog_halo), 1, out);
  }
  fclose(in);
  fclose(out);
  return 0;
}
