#ifndef _CORR_H_
#define _CORR_H_

#include <stdint.h>

struct reduced_halo {
  float pos[3];
  float sm;
  int nq;
};

struct point {
  float pos[3];
#ifndef MINIMAL_CORR_STRUCT
  int div; 
#endif /*!def MINIMAL_CORR_STRUCT*/
};

#define DIVISIONS 4
#define TOTAL_DIV (DIVISIONS*DIVISIONS)
#define MIN_DIST -3
#define MAX_DIST 1
#define DIST_BPDEX 5
#define NUM_BINS (((MAX_DIST-MIN_DIST) * DIST_BPDEX)+1)
#define INV_DIST_BPDEX (1.0/(double)DIST_BPDEX)

#define F_LOG2_DIVS 4
#define F_TOTAL_DIV (2<<(F_LOG2_DIVS))
#define F_INV_MIN_DIST2 (1e6)
#define F_MIN_DIST (sqrt(1.0/F_INV_MIN_DIST2))
#define F_NUM_BINS (27)
#define F_MAX_DIST (F_MIN_DIST*pow(2.0, (F_NUM_BINS+1.0)/2.0))
#define F_NUM_N2 (500)


#endif /* _CORR_H_ */
