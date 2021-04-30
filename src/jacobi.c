#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "all_smf.h"
//#define double long double

//#define NUM_PARAMS 5

void print_matrix(double matrix[][NUM_PARAMS]) {
  int i,j;
  for (i=0; i<NUM_PARAMS; i++) {
    for (j=0; j<NUM_PARAMS; j++) {
      printf("%+3e ", matrix[i][j]);
    }
    printf("\n");
  }
}

//Algorithm from the Wikipedia entry on Jacobi decomposition
void jacobi_decompose(double cov_matrix[][NUM_PARAMS], double *eigenvalues, double orth_matrix[][NUM_PARAMS]) {
  int i,j,k,l;
  int max_col[NUM_PARAMS];
  int changed[NUM_PARAMS]; //To keep track of eigenvalues changing
  int n = NUM_PARAMS;
  int max_row;
  int state = 0;
  double max, c, s, t, u, a;
  int count = 0;

  //Set up the maximum value cache
  for (i=0; i<n; i++) {
    max=0;
    max_col[i] = i+1;
    for (j=i+1; j<n; j++) {
      if (fabs(cov_matrix[i][j])>max) {
	max_col[i] = j;
	max = fabs(cov_matrix[i][j]);
      }
    }
  }

  //Set up the orthogonal matrix as the identity:
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      orth_matrix[i][j] = (i==j) ? 1 : 0;

  //Copy the diagonal values to the eigenvalue array:
  for (i=0; i<n; i++) {
    eigenvalues[i] = cov_matrix[i][i];
    changed[i] = 1;
    for (j=0; j<n; j++)
      if (j!=i && cov_matrix[i][j]) break;
    if (j==n) {
      state--;
      changed[i] = 0;
    }
  }

  //Sweep time: iterate until the eigenvalues stop changing.
  state += n;
  while (state) {
    count++;
    //Find the largest nonzero element in the matrix:
    max = fabs(cov_matrix[0][max_col[0]]);
    max_row = 0;
    for (i=1; i<n-1; i++) {
      if (fabs(cov_matrix[i][max_col[i]])>max) {
	max = fabs(cov_matrix[i][max_col[i]]);
	max_row = i;
      }
    }
    k = max_row; l = max_col[k];
    max = cov_matrix[k][l];
    if (max==0) break;

    //Calculate the Jacobi rotation matrix
    //Tan 2phi = 2S_kl / (S_kk - S_ll)
    a = (eigenvalues[l] - eigenvalues[k]) * 0.5;
    t = fabsl(a) + sqrtl(max*max + a*a);
    s = sqrtl(max*max + t*t);
    c = t/s;
    s = max/s;
    t = max*max / t;
    if (a<0) {
      s = -s; t = -t;
    }

    //Update eigenvalues
#define UPDATE(x,y)							\
    a = eigenvalues[x];							\
    eigenvalues[x] += y;						\
    if (changed[x] && (fabs(y) < 1e-6*fabs(eigenvalues[x]))) { /*Eigenvalue didn't change*/ \
      changed[x]=0;							\
      state--;								\
    } else if (!changed[x] && (fabs(y) > 1e-6*fabs(eigenvalues[x]))) { /*Egval did change*/ \
      changed[x] = 1;							\
      state++;								\
    }
      
    UPDATE(k, -t);
    UPDATE(l, t);

    //Update covariance matrix:
    cov_matrix[k][l] = 0;

#define ROTATE(m,w,x,y,z,r)  /*Perform a Jacobi rotation*/	\
      t = m[w][x]; u = m[y][z];					\
      m[w][x] = t*c - s*u;					\
      m[y][z] = s*t + c*u;				        \
      if (r) {							\
	if (fabs(m[w][x])>fabs(m[w][max_col[w]]))		\
	  max_col[w] = x;					\
	if (fabs(m[y][z])>fabs(m[y][max_col[y]]))		\
	  max_col[y] = z;					\
      }
    

    for (i=0; i<k; i++) { ROTATE(cov_matrix,i,k,i,l,1); }
    for (i=k+1; i<l; i++) { ROTATE(cov_matrix,k,i,i,l,1); }
    for (i=l+1; i<n; i++) { ROTATE(cov_matrix,k,i,l,i,1); }
    for (i=0; i<n; i++) { ROTATE(orth_matrix,k,i,l,i,0); }
    
#undef ROTATE
#undef UPDATE
  }

  // Restore the matrix to its original form:
  for (k=0; k<n-1; k++) {
    for (l=k+1; l<n; l++) {
      cov_matrix[k][l] = cov_matrix[l][k];
    }
  }
}


void matrix_multiply(double m[][NUM_PARAMS], double *in, double *out) {
  int i,j;
  for (i=0; i<NUM_PARAMS; i++) out[i] = 0;
  for (i=0; i<NUM_PARAMS; i++)
    for (j=0; j<NUM_PARAMS; j++)
      out[i] += m[i][j]*in[j];
}

void vector_madd(double *a, double b, double *c) {
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    a[i] += b*c[i];
}

int test(void) {
  /*  double S[NUM_PARAMS][NUM_PARAMS] = {{ 4,   -30,    60,   -35, 0},
				      {-30,  300,  -675,   420, 0},
				      { 60, -675,  1620, -1050, 0},
				      {-35,  420, -1050,   700, 0},
				      {  0,    0,     0,     0, 1}};
  */
  double S2[NUM_PARAMS][NUM_PARAMS] = {{ 5, -1, 0, 0, 0},
				       { -1, 2, 0, 0, 0},
				       { 0, 0, 1, 0, 0},
				       { 0, 0, 0, 1, 0},
				       { 0, 0, 0, 0, 1}};
  double O[NUM_PARAMS][NUM_PARAMS];
  double E[NUM_PARAMS];
  double IN[NUM_PARAMS], OUT[NUM_PARAMS];
  int i, j;
  jacobi_decompose(S2,E,O);
  for (i=0; i<NUM_PARAMS; i++) {
    printf("Eigenvalue[%d]: %f\n", i, E[i]);
    for (j=0; j<NUM_PARAMS; j++) {
      printf("Eigv[%d]: %f\n", j, O[i][j]);
      IN[j] = O[i][j];
    }
    matrix_multiply(S2,IN,OUT);
    for (j=0; j<NUM_PARAMS; j++) {
      printf("Eigv test[%d]: %f %f\n", j, OUT[j], O[i][j]*E[i]);
    }
    printf("\n");
  }
  return 0;
}
