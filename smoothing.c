#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>
#include "check_syscalls.h"


double cdf_diff(double x, double dx) {
  return (0.5*(erf(x+dx)-erf(x-dx)));
}

extern void cloud_in_cell(double val, double weight, double *bins, double min, int64_t bpunit, int64_t nbins);
extern void cloud_cell_indices(float val, double min, int64_t bpunit, int64_t nbins, int64_t *b, int64_t *b2, float *f, float *f2);

// void cloud_in_cell(double val, double weight, double *bins, double min, int64_t bpunit, int64_t nbins) {
// // void cloud_in_cell(double val, double weight, double *bins, double min, int64_t bpunit, int64_t nbins) {
//   float f = (val-min)*bpunit;
//   int64_t b = f-0.5; //Where left edge of cloud falls
//   int64_t b2 = b+1;
//   float f2 = f+0.5-b2; //Fraction of cloud in bin 2
//   f = 1.0-f2;
//   if (b>0 && b<nbins) bins[b] += f*weight;
//   if (b2>0 && b2<nbins) bins[b2] += f2*weight;
// }

// void cloud_cell_indices(float val, double min, int64_t bpunit, int64_t nbins, int64_t *b, int64_t *b2, float *f, float *f2) {
// // void cloud_cell_indices(float val, double min, int64_t bpunit, int64_t nbins, int64_t *b, int64_t *b2, float *f, float *f2) {
//   *f = (val-min)*bpunit;
//   *b = *f-0.5; //Where left edge of cloud falls
//   *b2 = *b+1;
//   *f2 = (*f)+0.5-*b2; //Fraction of cloud in bin 2
//   *f = 1.0-(*f2);
// }

void smooth_bins(double *input, double *output, double smoothing, double in_min, int64_t in_bpunit, int64_t in_nbins,
		 double out_min, int64_t out_bpunit, int64_t out_nbins, int64_t ssfr_correction, int64_t conserve_probability) {
  int64_t i, j;
  double *results = NULL;
  check_realloc_s(results, in_nbins, sizeof(double));
  if (!((in_bpunit % out_bpunit)==0)) {
    fprintf(stderr, "[Error] %"PRId64" %"PRId64"\n", in_bpunit, out_bpunit);
  }
  assert((in_bpunit % out_bpunit)==0);  // Make sure that each output bin contains an integer number of input bins

  // Make sure that the input is a superset of the output
  assert(in_min <= out_min);
  assert(in_min*in_bpunit+in_nbins >= out_min*in_bpunit+out_nbins*(in_bpunit/out_bpunit));

  //Make sure that input and output bins are aligned.
  double bin_difference = (out_min-in_min)*in_bpunit;
  int64_t bin_offset = bin_difference+1e-5;
  assert(fabs(bin_difference-bin_offset) < 1e-4);

  //Smooth if necessary
  if (smoothing>0) {
    double *stencil = NULL, *stencil_cdf = NULL;
    check_realloc_s(stencil, in_nbins*2+1, sizeof(double));
    memset(results, 0, sizeof(double)*in_nbins);
    double smoothing_const = 1.0/(sqrt(2.0)*smoothing*in_bpunit); //Because the CDF is erf(x/(width*sqrt(2)))
    double dx = 0.5*smoothing_const;
    for (i=0; i<in_nbins*2+1; i++) {
      double x = i-in_nbins;
      stencil[i] = cdf_diff(x*smoothing_const, dx);
    }
    if (ssfr_correction) {
      for (i=0; i<in_nbins*2+1; i++) {
	double x = (i-in_nbins)/(double)in_bpunit;
	stencil[i] *= pow(10,-0.35*x);  //Assume that obs. scatter in SFR \propto scatter in M_*^0.65
	                                //Makes it extra easy to deal with Moustakas quenched threshold.
      }
    }
    if (conserve_probability) {
      check_realloc_s(stencil_cdf, in_nbins*2+1, sizeof(double));
      stencil_cdf[0] = 0;
      for (i=1; i<in_nbins*2+1; i++) stencil_cdf[i] = stencil_cdf[i-1]+stencil[i-1];
    }
    for (i=0; i<in_nbins; i++) {
      if (!input[i]) continue;
      double val = input[i];
      double *ts = stencil+in_nbins-i;
      for (j=0; j<in_nbins; j++) results[j] += val*ts[j];
      if (conserve_probability) {
	results[0] += stencil_cdf[in_nbins-i];
	results[in_nbins-1] += 1.0-stencil_cdf[2*in_nbins-i-1];
      }
    }
    free(stencil);
    if (conserve_probability) free(stencil_cdf);
  } else {
    for (i=0; i<in_nbins; i++) results[i] = input[i];
  }

  //Now, rebin:
  memset(output, 0, sizeof(float)*out_nbins);
  for (i=0; i<out_nbins; i++) {
    j = i*(in_bpunit/out_bpunit) + bin_offset;
    int64_t max_j = j+(in_bpunit/out_bpunit);
    double v = 0;
    for (; j<max_j; j++) v += results[j];
    output[i] = v;
  }

  free(results);
}
