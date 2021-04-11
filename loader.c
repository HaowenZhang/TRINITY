#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>

#define A_UV(x) ((x).t_tdyn)
struct catalog_halo {
  int64_t id, descid, upid;
  int32_t flags; float uparent_dist;
  float pos[6], vmp, lvmp, mp, m, v, r;
  float rank1, rank2, ra, rarank, t_tdyn;
  float sm, icl, sfr, obs_sm, obs_sfr, obs_uv;
};

#define IGNORE_FLAG 16

/*
Field explanations:
**Note that halo masses are in Msun/h and stellar masses/SFRs are in Msun.
ID: Unique halo ID
DescID: ID of descendant halo (or -1 at z=0).
UPID: -1 for central halos, otherwise, ID of largest parent halo
Flags: Mostly internal UniverseMachine info.  However, halos with bit 4 set in the flags (i.e., flags & (2**4) is true) should be ignored.
Uparent_Dist: Ignore
pos[6]: (X,Y,Z,VX,VY,VZ)
X Y Z: halo position (comoving Mpc/h)
VX VY VZ: halo velocity (physical peculiar km/s)
M: Halo mass (Bryan & Norman 1998 virial mass, Msun/h)
V: Halo vmax (physical km/s)
MP: Halo peak historical mass (BN98 vir, Msun/h)
VMP: Halo vmax at the time when peak mass was reached.
R: Halo radius (BN98 vir, comoving kpc/h)
Rank1: halo rank in Delta_vmax (see UniverseMachine paper)
Rank2, RA, RARank: Ignore
A_UV: UV attenuation (mag)
SM: True stellar mass (Msun)
ICL: True intracluster stellar mass (Msun)
SFR: True star formation rate (Msun/yr)
Obs_SM: observed stellar mass, including random & systematic errors (Msun)
Obs_SFR: observed SFR, including random & systematic errors (Msun/yr)
Obs_UV: Observed UV Magnitude (M_1500 AB)
*/

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
  
  num_halos = file_stat.st_size /  sizeof(struct catalog_halo);
  halos = malloc(sizeof(struct catalog_halo)*num_halos);
  if (!halos) {
    fprintf(stderr, "[Error] Could not allocate memory for halos.");
    exit(EXIT_FAILURE);
  }
  
  if (fread(halos, sizeof(struct catalog_halo), num_halos, input)<0) {
    fprintf(stderr, "[Error] Could not read file %s!\n", filename);
    exit(EXIT_FAILURE);    
  }
  fclose(input);

  return halos;
}

int main(int argc, char **argv) {
  struct stat file_stat;
  FILE *input=NULL;
  char *filename = "sfr_catalog_0.055623.bin";
  struct catalog_halo *halos = NULL;
  int64_t i, num_halos=0;

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
  
  num_halos = file_stat.st_size /  sizeof(struct catalog_halo);
  halos = malloc(sizeof(struct catalog_halo)*num_halos);
  if (!halos) {
    fprintf(stderr, "[Error] Could not allocate memory for halos.");
    exit(EXIT_FAILURE);
  }
  
  if (fread(halos, sizeof(struct catalog_halo), num_halos, input)<0) {
    fprintf(stderr, "[Error] Could not read file %s!\n", filename);
    exit(EXIT_FAILURE);    
  }
  fclose(input);

  for (i=0; i<num_halos; i++) {
    if (halos[i].flags & IGNORE_FLAG) continue; /* Ignore halos with IGNORE_FLAG set */
    /* Rest of program logic here. */
  }
  return 0;
}
