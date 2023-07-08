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
  printf("#ID DescID X Y Z VX VY VZ M V MP VP R Sat Acc_Q Acc_QNE A0.5 A_4p A12 SM SFR OBS_SFR Acc_Scale Spin_Acc V_acc a_Vhalf V_mpeak C_acc SSFR\n");
  for (i=0; i<num_p; i++) {
    struct catalog_halo *tp = p+i;
    float ssfr = 0;
    if (tp->sm > 0) ssfr = tp->obs_sfr / tp->sm;
    printf("%"PRId64" %"PRId64" %f %f %f %f %f %f %.3e %f %.3e %f %f %.3e %.3e %.3e %f %f %f %e %e %e %f %f %f %f %f %f %e\n",
	   tp->id, tp->descid, tp->pos[0], tp->pos[1], tp->pos[2], tp->pos[3], tp->pos[4], tp->pos[5], tp->m, tp->v, tp->mp, tp->vp, tp->r, tp->sat, tp->acc_q, tp->acc_qne, tp->a_half, tp->a_4p, tp->a_12, tp->sm, tp->sfr, tp->obs_sfr, tp->acc_scale, tp->spin_acc, tp->va, tp->a_v05, tp->vmp, tp->c, ssfr);
  }
  munmap(p, new_length);
  return 0;
}
