#ifndef __UNIVERSE_TIME_H__
#define __UNIVERSE_TIME_H__

void init_time_table(double Om, double h0);
double scale_to_time(double scale);
double scale_to_years(double scale);
double years_to_scale(double years);

#endif
