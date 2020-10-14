#ifndef JACOBI_H
#define JACOBI_H

#include "all_smf.h"

void jacobi_decompose(double cov_matrix[][NUM_PARAMS], double *eigenvalues, double orth_matrix[][NUM_PARAMS]);

void vector_madd(double *a, double b, double *c); //Performs a += b*c

#endif
