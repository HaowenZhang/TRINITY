#include <math.h>
#include <stdio.h>
#include "universe_time.h"
#include "mah.h"

int main(void) {
	float Om=0.27, h=0.7, halo_mass = 13, scale = 1;
	init_time_table(Om, h);
	int i,j;
	for (i=0; i<10; i++) {
	  double z = i;
	  double a = 1.0/(i+1.0);
	  char buffer[1024];
	  snprintf(buffer, 1024, "results/halo_specific_accretion_z%d.dat", (int)z);
	  FILE *output = fopen(buffer, "w");
	  fprintf(output, "#M SpecificAcc\n");
	  for (j=90; j<131; j++) {
	    double m = j/10.0;
	    fprintf(output, "%e %e\n", pow(10, m), ma_rate_avg_mnow(m, a)/pow(10, m));
	  }
	  fclose(output);
	}
	return 0;
}
