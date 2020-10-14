#ifndef _SMOOTHING_H_
#define _SMOOTHING_H_

//Input / output bins must be aligned
//Input array must be superset of output
//BPunit is the # of bins per unit (e.g., dex, mag, etc.)
//Smoothing is the width of the normal scatter applied
void smooth_bins(double *input, double *output, double smoothing, double in_min, int64_t in_bpunit, int64_t in_nbins,
		 double out_min, int64_t out_bpunit, int64_t out_nbins, int64_t ssfr_correction, int64_t conserve_probability);

// void cloud_in_cell(double val, double weight, double *bins, double min, int64_t bpunit, int64_t nbins);
// void cloud_cell_indices(float val, double min, int64_t bpunit, int64_t nbins, int64_t *b, int64_t *b2, float *f, float *f2);


inline void cloud_in_cell(double val, double weight, double *bins, double min, int64_t bpunit, int64_t nbins) {
// void cloud_in_cell(double val, double weight, double *bins, double min, int64_t bpunit, int64_t nbins) {
  float f = (val-min)*bpunit;
  int64_t b = f-0.5; //Where left edge of cloud falls
  int64_t b2 = b+1;
  float f2 = f+0.5-b2; //Fraction of cloud in bin 2
  f = 1.0-f2;
  if (b>0 && b<nbins) bins[b] += f*weight;
  if (b2>0 && b2<nbins) bins[b2] += f2*weight;
}

inline void cloud_cell_indices(float val, double min, int64_t bpunit, int64_t nbins, int64_t *b, int64_t *b2, float *f, float *f2) {
// void cloud_cell_indices(float val, double min, int64_t bpunit, int64_t nbins, int64_t *b, int64_t *b2, float *f, float *f2) {
  *f = (val-min)*bpunit;
  *b = *f-0.5; //Where left edge of cloud falls
  *b2 = *b+1;
  *f2 = (*f)+0.5-*b2; //Fraction of cloud in bin 2
  *f = 1.0-(*f2);
}


double cdf_diff(double x, double dx);

#endif /* _SMOOTHING_H_ */
