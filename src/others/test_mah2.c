#include <math.h>
#include <stdio.h>
#include "universe_time.h"
#include "mah.h"

int main(void) {
	float Om=0.27, h=0.7, halo_mass = 13, scale = 1;
	init_time_table(Om, h);
	int i;
	for (i=0; i<15; i++) {
	  double a = 1.0/(i+1.0);
	  //printf("Mass history for 10^13 Msun halo at a=%g: %e Msun; %e Msun\n", a, pow(10,m_evolution_avg(13, 1, a)), pow(10, m_at_a(13, a)));
	  printf("%d %f\n", i, m_evolution_avg(14, 1, a));
	}
	return 0;
}
