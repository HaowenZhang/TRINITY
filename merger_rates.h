#ifndef _MERGER_RATES_H_
#define _MERGER_RATES_H_

//Log units of mass, ratio
//Returns merger rate / halo / unit redshift / unit log ratio
//Mass in units of Msun (no h).
double merger_rate_macc(double mass, double ratio, double scale);
double merger_rate_mpeak(double mass, double ratio, double scale);

#endif /* _MERGER_RATES_H_ */
