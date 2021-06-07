#ifndef _MAKE_GRAPHS_H_
#define _MAKE_GRAPHS_H_
#include <inttypes.h>
#include "all_smf.h"

void print_sfh(int64_t count, struct smf_fit input);
int64_t find_step(double scale);
static void read_params(char *buffer, double *data, int max_n);
struct smf_fit parse_smf_fit(char *buffer);

#endif /* _MAKE_GRAPHS_H_ */
