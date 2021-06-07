#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "observations.h"
#include "mf_cache.h"
#include "integrate.h"
#include "integrate2.h"
#include "distance.h"
#include "smf.h"
#include "expcache2.h"
#include "universe_time.h"
#include "integrate2.c"
#include "all_smf.h"
#include "calc_sfh.h"
#include "mlist.h"
#include "check_syscalls.h"
#include <omp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

// The parameters to calculate the correction for Compton-thick AGNs.
#define BETA 0.24
#define psi_min 0.2
#define psi_max 0.84
#define psi43750 0.43
#define a1 0.48

#define LOG10_4 0.6020599913279624

extern int no_z_scaling;
extern float z_max;
extern float z_min;



#define exp10fc(x) doexp10(x)

// The function to do 10^x.
static double doexp10(double x) 
{
  double a = exp(M_LN10*x);
  return a;
}

// The functions to calculate the fraction of Compton thick 
// obscured AGNs. See Ueda et al. (2014).
double psi4375(double z)
{
  return z < 2? psi43750 * pow((1 + z), a1) : psi43750 * pow((1 + 2), a1);
}

double psi(double Lx, double z)
{
  double max1 = psi4375(z) - BETA * (Lx - 43.75);
  max1 = max1 > psi_min ? max1 : psi_min;
  return psi_max > max1 ? max1 : psi_max;
}

// Find the X-ray luminosity given the bolometric luminosity.
// All the luminosities are in log10 units. See Ueda et al. (2014.)
double find_Lx_at_Lbol(double Lbol)
{
  double BC = 10.83 * pow(exp10(Lbol) / (1e10 * 3.839e33), 0.28) +
                            6.08 * pow(exp10(Lbol) / (1e10 * 3.839e33), -0.02);
  return Lbol - log10(BC);
}

// Find the bolometric luminosity given the X-ray luminosity.
// All the luminosities are in log10 units. See Ueda et al. (2014).
double find_Lbol_at_Lx(double Lx)
{
  double dlx = Lx - 4.21484284e+01;
  double logBC = 8.82970680e-01 + log10(exp10(-3.31381964e-03 * dlx) + exp10(3.87414950e-01 * dlx));
  return Lx + logBC;
}

// The fraction of obscured AGNs, as a function of X-ray
// luminosity. This is used to correct for the active
// black hole mass functions from Kelly & Shen (2013).
double F_obs(double Lx)
{
  if (Lx < 43) return 0.985;
  return 0.56 + 1.0 / M_PI * atan((43.89 - Lx) / 0.46);
}

int use_obs_psf = 1; // The flag indicating if the scatter in observed stellar masses
                    // is accounted for in the calculation of observables.
struct z_mf_cache *all_mf_caches = NULL; // The cache of halo mass functions.
double gauss_cache[GAUSS_CACHE_SIZE]; // The cache for the standard normal distribution.

// Initialize the cache of standard normal distributions.
void init_gauss_cache(void) 
{
  int i;
  double m;
  for (i=0; i<GAUSS_CACHE_SIZE; i++) 
  {
    m = (double)i*GAUSS_CACHE_STEP + GAUSS_CACHE_MIN;
    gauss_cache[i] = exp(-0.5*m*m)/sqrt(2.0*M_PI);
  }
}

// Calculate the normalization of a Schechter function with a power-law index
// alpha, within the interval [lower, upper]. Note that the input lower/upper
// are in log10 units.
// Note tht this function is effectively deprecated, since we adopt a double
// power-law Eddington ratio distribution.
double schechter_norm(double alpha, double lower, double upper, struct smf_fit *fit) 
{
  lower = exp10fc(lower);
  upper = exp10fc(upper);
  gsl_sf_result rl, rh;

  // The integral can be obtained by taking the difference of two incomplete
  // gamma functions.
  if (gsl_sf_gamma_inc_e(alpha, lower, &rl) || gsl_sf_gamma_inc_e(alpha, upper, &rh)) 
  {
    //error happened...
    if (fit) 
    {
      INVALIDATE(fit, "Invalid arguments to incomplete gamma function.");
    }
    return 1;
  }
  return((rl.val-rh.val)/log(10.0));
}

// Calculate the normalization of a double power-law: 1 / (x**a + x**b).
double doublePL_norm(double a, double b, double lower, double upper, struct smf_fit *fit) 
{ 
  // Ensure that a is bigger than b.
  if (a < b)
  {
    double tmp = a;
    a = b;
    b = tmp;
  }

  // Important: P(bher) is per log(eta), so to convert it to an integral with deta,
  // we need to add 1 to both of the power indices.
  a += 1;
  b += 1;

  // Generally, the definite integral of this double power-law is a confluent 
  // hypergeometric function. I don't wanna get into the weeds, and so implemented
  // a numerical integral below.
  int n_pts = 200; // # of points in the numerical integral
  double scale = 1.0 / n_pts * (upper - lower); // scale is the multiplicative factor
                                                // connecting the midway point of
                                                // consecutive small intervals.
  scale = exp10fc(scale);
  double x1 = exp10fc(lower);
  double x2 = exp10fc(upper);
  double val = 0;
  double x_left = x1;
  double x_right = x_left * scale;
  double x_center = 0.5 * (x_left + x_right);
  double dx = x_right - x_left;
  double x_center_a = pow(x_center, a);
  double x_center_b = pow(x_center, b);
  double scale_a = pow(scale, a); // Precalculate the scale factor of the two power-laws,
  double scale_b = pow(scale, b); // so that we don't have to invoke pow() every time.

  for (int i = 0; i < n_pts; i++)
  {    
    val += 1 / (x_center_a + x_center_b) * dx;
    dx *= scale;
    x_center_a *= scale_a;
    x_center_b *= scale_b;
  }
  
  // Also we need a 1/(ln10) in the normalization.
  return val / log(10.0);
}


// Calculate the inverse of the average Eddington ratio, 
// given the Schechter power-law index, within the
// interval [lower, upper].
double schechter_inv_avg(double alpha, double lower, double upper, struct smf_fit *fit) 
{
  double norm = schechter_norm(alpha, lower, upper, fit);
  double norm2 = schechter_norm(alpha+1, lower, upper, fit);
  return norm/norm2;
}

// Calculate the fraction of schechter distribution 
// (power-law index alpha, normalization norm)
// that is above a certain threshold (thresh).
double schechter_frac_above_thresh(double thresh, double alpha, double norm, struct smf_fit *fit) 
{
  if (thresh < BHER_EFF_MIN) thresh = BHER_EFF_MIN;
  if (thresh >= BHER_EFF_MAX) return 0;
  return (schechter_norm(alpha, thresh, BHER_EFF_MAX, fit)*norm);
}

// Calculate the fraction of double power-law
// distribution (power-law indices alpha and delta,
// normalization norm) that is above a certain,
// threshold (thresh).
double doublePL_frac_above_thresh(double thresh, double alpha, double delta, double norm, struct smf_fit *fit) 
{
  if (thresh < BHER_EFF_MIN) thresh = BHER_EFF_MIN;
  if (thresh >= BHER_EFF_MAX) return 0;
  return (doublePL_norm(alpha, delta, thresh, BHER_EFF_MAX, fit)*norm);
}

// Calculate the Gaussian distribution value based on the 
// difference with the central value, scaled by the Gaussian
// spread (dm).
double syst_cache(double dm) 
{
  return (exp(-0.5*dm*dm)*(0.5*M_2_SQRTPI*M_SQRT1_2));
}

// Calculate the inverse of the total scatter,
// which is a quadratic sum of the scatter in
// observed stellar masses around the true masses
// (obs_scatter), and the intrinsic scatter around,
// e.g., the stellar mass--halo mass relation.
double gen_inv_sigma(double obs_scatter, double scatter) 
{
  if (!scatter && (!use_obs_psf || !obs_scatter)) return 0;
  double gauss_sigma = sqrt(obs_scatter*obs_scatter + scatter*scatter);
  return 1.0/gauss_sigma;
}

// Evaluate the point spread function at a certain
// mass difference (delta_m), given the inverse
// of total scatter (gauss_inv_sigma).
double evaluate_psf(double delta_m, double gauss_inv_sigma) 
{
  if (!gauss_inv_sigma) return (delta_m ? 0 : 1);
  return(gauss_inv_sigma*syst_cache(delta_m*gauss_inv_sigma));
}

// Interpolation of list over stellar mass.
// This function is used only when list contains
// stellar masses.
double _interp_from_sm(double sm, double *list, int64_t n, char *ok) 
{
  double f = (sm-SM_MIN)*((double)SM_BPDEX)+SM_EXTRA;
  int64_t b = f;
  if (b >= SM_EXTRA+SM_BINS-1) 
    {
      // printf("the sm=%f value exceeded the bin number. \n", sm);
      return 1e-17;
    }
  if (f < 0) return 1e-17;
  f-=b;
  if (n>-1 && !(ok[b] && ok[b+1])) return calc_smf_at_sm(n, sm);
  return (list[b] + f*(list[b+1]-list[b]));
}

// Interpolation of list over stellar mass.
// This function is used only when the list
// contains SFRs.
double _interp_from_sfr(double sm, double *list, int64_t n) 
{
  double f = (sm-SM_MIN)*((double)SM_BPDEX)+SM_EXTRA;
  int64_t b = f;
  if (b >= SM_EXTRA+SM_BINS-1) 
    {
      // printf("the sm=%f value exceeded the bin number. \n", sm);
      return 1e-17;
    }
  if (f < 0) return 1e-17;
  f-=b;
  // if (n>-1 && !(ok[b] && ok[b+1])) return calc_smf_at_sm(n, sm);
  return (list[b] + f*(list[b+1]-list[b]));
}

// Interpolation of list over stellar mass.
// This function is used only when the list
// contains galaxy quenched fractions.
double _interp_from_sfrac(double sm, double *list, int64_t n) 
{
  double f = (sm-SM_MIN)*((double)SM_BPDEX)+SM_EXTRA;
  int64_t b = f;
  if (b >= SM_EXTRA+SM_BINS-1) 
    {
      // printf("the sm=%f value exceeded the bin number. \n", sm);
      return 1e-17;
    }
  if (f < 0) return 1e-17;
  f-=b;
  // if (n>-1 && !(ok[b] && ok[b+1])) return calc_smf_at_sm(n, sm);
  return (list[b] + f*(list[b+1]-list[b]));
}

// Interpolation of list over stellar mass.
// This function is used only when the list
// contains UV magnitudes.
double _interp_from_uv(double uv, double *list, int64_t n, char *ok) 
{
  double f = (uv-UV_MIN)*((double)UV_BPMAG)+UV_EXTRA;
  int64_t b = f;
  if (b >= UV_EXTRA+UV_BINS-1) 
    {
      // printf("the sm=%f value exceeded the bin number. \n", sm);
      return 1e-17;
    }
  if (f < 0) return 1e-17;
  f-=b;
  if (n>-1 && !(ok[b] && ok[b+1])) 
  {
  	double result = calc_uvlf_at_uv(n, uv);
	return result;
  }
  return (list[b] + f*(list[b+1]-list[b]));
}

// Interpolation of list over stellar mass.
// This function is used only when the list
// contains the Gaussian spread of UV
// magnitudes.
double _interp_from_std_uv(double uv, double *list, int64_t n) 
{
  double f = (uv-UV_MIN)*((double)UV_BPMAG)+UV_EXTRA;
  int64_t b = f;
  if (b >= UV_EXTRA+UV_BINS-1)
    {
	      return 1e-17;
    }
  if (f < 0) return 1e-17;
  f-=b;
 
  return (list[b] + f*(list[b+1]-list[b]));
}
  

// Calculate the ***observed*** bulge mass
// given the ***observed*** stellar mass
// and scale factor (a= 1 / (1 + z)).
double bulge_mass(double sm, double a) 
{
  double z_mul = 1.0 - 0.5*(1.0-a);
  return sm+log10(z_mul/(1.0+exp(-1.12882*(sm-10.1993))));
}

// Calculate the contribution of a certain stellar mass
// (sm) to the stellar mass function (SMF) at a pre-
// defined stellar mass, extra_data->sm.
double step_integral_helper(double sm, void *extra_data) 
{
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;

  // Evaluate the point spread function, which specifies
  // how much contribution from sm is scattered to ih->sm.
  double psf1 = evaluate_psf(delta_m, ih->gauss_inv_sigma);


  // Apply further corrections if needed.
  if (ih->corr)
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
  if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
  if (!isfinite(val)) 
  {
    val = 1e-17;
  }
  return (psf1*val);
}

// The new step_integral_helper for ssfr. This is because we have to separate the
// lists to be interpolated (i.e., SFR*SMF) for SF and Q galaxies.
double step_integral_helper_ssfr(double sm, void *extra_data) 
{
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;

  // We need two PSFs in this function, because both the contributions from star-forming
  // and quenched galaxies should be accounted for.
  double sfrac = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  double psf1 = (1 - sfrac) * evaluate_psf(delta_m, ih->gauss_inv_sigma);
  double psf2 = psf1 / (1 - sfrac) * sfrac;

  if (ih->corr)
  {
    // psf *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf2 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  }
  // Get the average SFRs for both star-forming and quenched galaxies.
  double val1 = _interp_from_sfr(sm, ih->sfr_sm_q, ih->n);
  double val2 = _interp_from_sfr(sm, ih->sfr_sm_sf, ih->n);
  if (ih->f > -1) val1 += ih->f*(_interp_from_sfr(sm, ih->sfr_sm_q2, ih->n2)-val1);
  if (ih->f > -1) val2 += ih->f*(_interp_from_sfr(sm, ih->sfr_sm_sf2, ih->n2)-val2);
  if ((!isfinite(val1)) || (!isfinite(val2))) 
  {
    val1 = val2 = 1e-17;
  }
  // Add up their contributions.
  return (psf1*val1 + psf2*val2);
}

// Calculate the contribution of a certain stellar mass
// (uv) to the UV luminosity function (UVLF) at a pre-
// defined UV magnitude, extra_data->sm (not extra_data->uv).
double step_integral_helper_uv(double uv, void *extra_data) 
{
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_uv = uv - ih->sm;

  // Here we use ih->sfrac_sm to contain the std_uvlf grids.
  double std_uv = _interp_from_std_uv(uv, ih->sfrac_sm, ih->n);
  double gauss_inv_sigma = gen_inv_sigma(0, std_uv);

  double psf = evaluate_psf(delta_uv, gauss_inv_sigma);


  double val = _interp_from_uv(uv, ih->list, ih->n, ih->ok);
  double val1 = val;
  double val2;
  if (ih->f > -1) val += ih->f*(_interp_from_uv(uv, ih->list2, ih->n2, ih->ok2)-val);
  val2 = _interp_from_uv(uv, ih->list2, ih->n2, ih->ok2);
  if (!isfinite(val)) 
  {
    val = 1e-17;
  }
  return (psf*val);
}

// Calculate the contribution of a certain stellar mass
// (sm) to the stellar mass function (SMF) at a pre-
// defined stellar mass, extra_data->sm. NOTE: This function
// only accounts for the quenched (Q) galaxies. This is gonna
// be used when calculating galaxy quenched fractions.
double step_integral_helper_Q(double sm, void *extra_data) 
{
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;
  double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);

  double sfrac1 = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  double psf1 = (1 - sfrac1) * evaluate_psf(delta_m, ih->gauss_inv_sigma);

  if (ih->corr)
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
  if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
  if (!isfinite(val)) 
  {
    val = 1e-17;
  }
  return (psf1*val);
}

// Calculate the contribution of a certain stellar mass
// (sm) to the stellar mass function (SMF) at a pre-
// defined stellar mass, extra_data->sm. NOTE: This function
// only accounts for the star-forming (SF) galaxies. 
// This is gonna be used when calculating galaxy quenched fractions.
double step_integral_helper_SF(double sm, void *extra_data) 
{
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;
  double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);

  double sfrac1 = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  double psf1 = sfrac1 * evaluate_psf(delta_m, ih->gauss_inv_sigma);

  if (ih->corr)
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
  if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
  if (!isfinite(val)) 
  {
    val = 1e-17;
  }
  return (psf1*val);
}

// Evaluate the observed stellar mass function
// or average observed specific SFR (SSFR) at 
// redshift z, ***observed*** stellar mass m.
double evaluate_from_step(double z, double sm, int smf) 
{
  // phi_corr is the correction in the normalizations of
  // either stellar mass function or average SSFR.
  double phi_corr=1, sm_true, sm_max;
  struct step_integral_helper_data ih;

  // Find out the corresponding snapshots for redshift z.
  calc_step_at_z(z, &(ih.step), &(ih.f));

  // Get the model parameters.
  struct smf smhm = steps[ih.step].smhm;

  // The correction between the average and median values due to
  // log-normal distribution.
  double scatter_corr = smhm.scatter*smhm.scatter*(0.5*M_LN10);

  // Correct for the systematic offset in stellar mass. See
  // Section 2.3 of Zhang et al. (2021).
  sm_true = (sm - smhm.mu);

  // Complain if the input stellar mass is unphysical.
  if (!(sm_true >= 0 && sm_true < 20)) 
  {
    fprintf(stderr, "Gah! %f %f %f %f %f %f %f %f %f %f %f\n", sm_true, sm, smhm.sm_0, smhm.v_1, smhm.delta, smhm.beta, smhm.gamma, smhm.lambda, smhm.scatter, smhm.mu, smhm.kappa);
  }
  
  // smf == 0 means that we're calculating the average SSFR. 
  // In this case we should apply the offset in SFR as
  // described in Section 2.3 of Zhang et al. (2021).
  if (!smf) 
  {
    phi_corr *= pow(10,-sm_true+smhm.kappa * exp(-0.5 * (z - 2) *(z - 2)))*smhm.ssfr_corr;
    // Also have to account for the correlation between stellar mass
    // and SFR at fixed halo mass.
    ih.corr = smhm.sfr_sm_corr;
  }

  // Otherwise we're calculating SMFs. No correction for 
  // the normalization should be applied.
  else
  {
    ih.corr = 0;
  }

  // Write the information into the help structure, step_integral_helper_data ih.
  if (!use_obs_psf) smhm.obs_scatter = 0;
  ih.gauss_inv_sigma = gen_inv_sigma(smhm.obs_scatter, smhm.scatter);
  ih.s_corr = scatter_corr;
  ih.sm = sm_true;
  ih.scatter = smhm.scatter;
  ih.obs_scatter = smhm.obs_scatter;
  ih.n = ih.n2 = -1;
  if (smf) 
  {
    // ih.list is the array we're going to interpolate on.
    ih.list = steps[ih.step].smf;
    ih.sfrac_sm = steps[ih.step].sfrac_sm;
    // ih.n is the step where we're doing the interpolation.
    ih.n = ih.step;
    // ih.ok shows if ih.list varies smoothly. See calc_smf_and_ssfr()
    // in calc_sfh.c.
    ih.ok = steps[ih.step].smf_ok;

    // If we're not looking at the last snapshot,
    // we also need to interpolate in the redshift
    // dimension.
    if (ih.f > -1) 
    {
      // Similar stuff but for the next snapshot
      ih.list2 = steps[ih.step+1].smf;
      ih.n2 = ih.step+1;
      ih.ok2 = steps[ih.step+1].smf_ok;
    }
  } 
  else // If we are calculating SSFRs
  {
    // ih.list is then the average SFR as
    // a function of stellar mass.
    ih.list = steps[ih.step].sfr_sm;
    // We also need star-forming fraction,
    // and SFRs for star-forming and quenched
    // galaxies.
    ih.sfrac_sm = steps[ih.step].sfrac_sm;
    ih.sfr_sm_sf = steps[ih.step].sfr_sm_sf;
    ih.sfr_sm_q = steps[ih.step].sfr_sm_q;

    // Same stuff for the next snapshot.
    if (ih.f > -1) 
    {
      ih.list2 = steps[ih.step+1].sfr_sm;
      ih.sfr_sm_sf2 = steps[ih.step+1].sfr_sm_sf;
      ih.sfr_sm_q2 = steps[ih.step+1].sfr_sm_q;
    }
  }

  // Due to various systematic offset and scatters,
  // the observed stellar mass function at any 
  // given stellar mass is a convolution of the
  // contribution from a broad range of intrinsic
  // stellar masses. And sm_max determines the 
  // maximum stellar mass whose contribution 
  // will be included when calculating the stellar 
  // mass functions.
  sm_max = sm_true+GAUSS_CACHE_MAX/ih.gauss_inv_sigma;
  if (sm_max > steps[ih.step].smhm.sm_max)
    sm_max = steps[ih.step].smhm.sm_max;

  double result;

  // No-scatter case.
  if (!smhm.scatter && (!use_obs_psf || !smhm.obs_scatter))
  {
    if (smf == 2) //smf == 2 means calculating quenched fractions
      {
        double smf_q = step_integral_helper_Q(sm_true - smhm.scatter, (void *)&ih);
        double smf_sf = step_integral_helper_SF(sm_true - smhm.scatter, (void *)&ih);
        result = smf_q / (smf_q + smf_sf);
        fprintf(stderr, "Interested SM: %.6f, sm_true: %.6f; smf_q: %.6e, smf_sf: %.6e\n",
              sm, sm_true, smf_q, smf_sf);
      }
    else
      result = (phi_corr*step_integral_helper(sm_true - smhm.scatter, (void *)&ih));
  }

  // If there are scatters in the model
  else 
  {
    double precision = PHI_INTEGRAL_PRECISION;
    // Require higher precision for very low redshifts. This is because
    // These data tend to have much smaller error bars.
    if (z <= 0.2) precision*=0.001;

    if (smf == 2) // Quenched fraction
    {

      double smf_q = (phi_corr * adaptiveGauss(step_integral_helper_Q, (void *)&ih,
          sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
            sm_max,
               precision,omp_get_thread_num()));
      double smf_sf = (phi_corr * adaptiveGauss(step_integral_helper_SF, (void *)&ih,
          sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
            sm_max,
               precision,omp_get_thread_num()));

      result = smf_q / (smf_q + smf_sf);
    }
    else if (smf == 1) // Stellar mass function
    {
      result = (phi_corr * adaptiveGauss(step_integral_helper, (void *)&ih,
          sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
            sm_max,
               precision,omp_get_thread_num())); 
    }

    else 
    {//smf == 0, i.e., SSFR calculation
      result = (phi_corr * adaptiveGauss(step_integral_helper_ssfr, (void *)&ih,
            sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
              sm_max,
                 precision,omp_get_thread_num()));
      double smf = evaluate_from_step(z, sm, 1);
      if (smf > 0) return result/smf;
      else return 0;
    }
  }
  return result;
}

// Calculate the observed UV luminosity function
// at a given observed UV magnitude uv and redshift
// z.
double evaluate_from_step_uv(double z, double uv) 
{
  double phi_corr=1;
  struct step_integral_helper_data ih;
  calc_step_at_z(z, &(ih.step), &(ih.f));
  struct smf smhm = steps[ih.step].smhm;

  // Similar info as in evaluate_from_step()
  ih.sm = uv;
  ih.n = ih.n2 = -1;
  ih.list = steps[ih.step].uvlf;
  ih.sfrac_sm = steps[ih.step].std_uvlf;
  ih.n = ih.step;
  ih.ok = steps[ih.step].uvlf_ok;
  
  // Same stuff for the next snapshot, if there
  // is a next snapshot.
  if (ih.f > -1) 
  {
    ih.list2 = steps[ih.step+1].uvlf;
    ih.n2 = ih.step+1;
    ih.ok2 = steps[ih.step+1].uvlf_ok;
  }
  
  // Determine the lower/upper limits for the integral.
  double uv_max = uv+UV_MAG_OFFSET;
  if (uv_max > steps[ih.step].smhm.uv_max)
    uv_max = steps[ih.step].smhm.uv_max;
  double uv_min = uv - UV_MAG_OFFSET;
  if (uv_min < steps[ih.step].smhm.uv_min)
    uv_min = steps[ih.step].smhm.uv_min;
 
  double result;
  double precision = PHI_INTEGRAL_PRECISION;
  // Higher precision for very low redshifts.
  // In fact this is very unlikely to be
  // triggered since UVLF at such low-z is 
  // very very difficult to calculate...
  if (z <= 0.2) precision*=0.001;

  result = (phi_corr * adaptiveGauss(step_integral_helper_uv, (void *)&ih,
      uv_min,
        uv_max,
           precision,omp_get_thread_num())); 
      
  return result;
}

// Calculate the average observed specific
// SFR (SSFR) as a function of ***observed***
// stellar mass (sm) and redshift (z).
double calc_ssfr(double sm, double z) 
{
  double ssfr = evaluate_from_step(z, sm, 0);
  if (ssfr < 1e-15) ssfr = 1e-15;
  return ssfr;
}

// Calculate the total (active + dormant)
// black hole mass function at a given black
// hole mass (m) and redshift (z). Note that
// this is an older implementation. With
// total BHMF already calculated at calc_bh_lum_distribution_full()
// in calc_sfh.c for each snapshot, it would
// be much faster to simply interpolate BHMF 
// over two consecutive snapshots. However,
// this one is kept here since it is used
// in MCMC, where the speed really matters.
double calc_bhmf(double m, double z) 
{
  int64_t step;
  double f;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  
  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-19;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;

  // Get the median BH mass and number densities of each halo mass
  // bin, by interpolating between the two consecutive snapshots.
  double bhm[M_BINS]={0}, t[M_BINS]={0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  int64_t s;
  // For higher precision, we further divide each halo mass bin
  // into 5 smaller bins.
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
  }

  nd /= 5.0; // Account for the fact that we split each halo
             // mass bin in 5.
  if (nd > 1e-19) return nd;
  return 1e-19;
}

double calc_bhmf_mh(double m, double z, double mh_low, double mh_high) 
{
  int64_t step;
  double f;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  
  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-19;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;

  // Get the median BH mass and number densities of each halo mass
  // bin, by interpolating between the two consecutive snapshots.
  double bhm[M_BINS]={0}, t[M_BINS]={0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  int64_t s;
  // For higher precision, we further divide each halo mass bin
  // into 5 smaller bins.

  int64_t mh_bpdex = 20;
  double inv_mh_bpdex = 1.0 / mh_bpdex;
  int64_t mh_bins = (mh_high - mh_low) * mh_bpdex;
  if (mh_bins < 1)
  {
    fprintf(stderr, "[mh_low, mh_high] interval is too small. The minimum difference should be %.6f dex.\n", inv_mh_bpdex);
    return 1e-19;
  }

  for (s=0; s<=mh_bins; s++) 
  {
    double mh = mh_low + (s + 0.5) * inv_mh_bpdex;
    f = (mh - M_MIN - 0.5 * INV_BPDEX) * BPDEX;
    i = f; f -= i;
    // printf("mh=%f, i=%d, f=%f\n", mh, i, f);
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
  }

  nd *= BPDEX * 1.0 / mh_bpdex; // Account for the fact that we split each halo mass bin into (mh_bpdex / BPDEX) sub-bins.
  if (nd > 1e-19) return nd;
  return 1e-19;
}

// Calculate the type I quasar mass function 
// at a given black hole mass (m) and redshift (z).
// See also: Kelly & Shen (2013)
double calc_bhmf_typeI(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }

  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  double bhm_f = (m - steps[step].bh_mass_min) * mbh_bpdex;
  int64_t bhm_b = bhm_f;
  bhm_f -= bhm_b;

  if (bhm_b >= MBH_BINS - 1) {bhm_b = MBH_BINS - 2; bhm_f = 1;}

  // We start with the total BHMFs that are pre-calculated.
  double bhmf_tot1 = steps[step].bhmf[bhm_b] + bhm_f * (steps[step].bhmf[bhm_b+1] - steps[step].bhmf[bhm_b]);
  double bhmf_tot2 = steps[step+1].bhmf[bhm_b] + bhm_f * (steps[step+1].bhmf[bhm_b+1] - steps[step+1].bhmf[bhm_b]);
  double bhmf_tot = bhmf_tot1 + f * (bhmf_tot2 - bhmf_tot1);

  // mass- and redshift-dependent AGN duty cycles.
  double f_mass = exp((m - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
  f_mass = f_mass / (1 + f_mass);
  double dc = steps[step].smhm.bh_duty * f_mass;
  if (dc < 1e-4) dc = 1e-4;

  // Obscured AGNs won't look like Type I, so we need to count 
  // how many of them are there and subtract them from
  // the total BHMFs.
  double corr_obs = 0;
  double tw = 0;
  for (i=0; i<LBOL_BINS; i++)
  {
    double lbol = LBOL_MIN + (i + 0.5) * LBOL_INV_BPDEX;
    double lx = find_Lx_at_Lbol(lbol);
    double f_obs = F_obs(lx); // The obscured fraction is a function of X-ray luminosity.
    corr_obs += (1 - f_obs) * steps[step].lum_dist_full[bhm_b*LBOL_BINS+i];
    tw += steps[step].lum_dist_full[bhm_b*LBOL_BINS+i];
  }  
  if (tw > 0) corr_obs /= tw; // The final Type I fraction.
  return bhmf_tot * corr_obs * dc;

}

// Calculate active black hole mass functions
// at a given black hole mass (m) and redshift (z).
// See also: Schulze & Wisotzki (2010), Schulze
// et al. (2015)
double calc_active_bhmf(double m, double z) 
{
  int64_t step;
  double f;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }

  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double s1 = (1.0-f)*steps[step].smhm.scatter*steps[step].smhm.bh_gamma + f*(steps[step+1].smhm.scatter)*steps[step+1].smhm.bh_gamma;
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double abhmf_shift = (1.0-f)*steps[step].smhm.abhmf_shift + f*(steps[step+1].smhm.abhmf_shift);
  if (z > 0.5) abhmf_shift = 0;
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;

  // Get the median BH mass and number densities of each halo mass
  // bin, by interpolating between the two consecutive snapshots.
  // Aside from these, we also have to get the active fraction (f_active)
  // that is calculated in calc_active_bh_fraction(), calc_sfh.c
  double bhm[M_BINS]={0}, t[M_BINS]={0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    double nactive = (1.0-f)*steps[step].f_active[i] + f*steps[step+1].f_active[i];
    t[i] = ndm*nactive;
  }

  // For higher precision, we further divide each halo mass bin
  // into 5 smaller bins.
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double dm = tbhm - (m + abhmf_shift);
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
  }
  nd /= 5.0; // Account for the fact that we split each halo
             // mass bin in 5.
  if (nd > 1e-15) return nd;
  return 1e-15;
}

// calculate the average BHAR as a function of 
// black hole mass (m) and redshift (z).
double calc_bhar_mbh(double m, double z) 
{
  int64_t step;
  double f;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }

  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0};
  double bhm_avg[M_BINS] = {0};
  double bhar[M_BINS] = {0};

  // Get the average, median BH masses, BH accretion rates,
  // and number densities of each halo mass bin, by interpolating 
  // between the two consecutive snapshots.
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // For higher precision, we further divide each halo mass bin
  // into 5 smaller bins.
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    
    // Average BH mass is used to account for the correlation between
    // BH mass and BH accretion rate. This correlation is induced by
    // our assumption that BHs of different masses at fixed ***halo mass***
    // share the same Eddington ratio distribution.
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    // With this assumption of Eddington ratio distribution, BH accretion
    // rates are enhanced by the same amount as the offset of BH mass
    // relative to the average mass.
    tbhar *= exp10fc(m - tbhm_avg); 
    double dm = tbhm - m;

    // The weight is determined by a log-normal distribution.
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
  }
  bhar_nd /= nd;
  if (bhar_nd > 1e-15) return bhar_nd;
  return 1e-15;
}

// calculate the average BHAR as a function of 
// stellar mass (Mstar) and redshift (z).
double calc_bhar_mstar(double m, double z) 
{

  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step].smhm.scatter);
  double mu = (1.0-f)*(steps[step].smhm.mu) + f*(steps[step].smhm.mu);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // Since we're looking at fixed stellar mass, the total scatter
  // is the one around the median stellar mass--halo mass relation.
  double s3 = s1*s1;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;

  // This function is very similar to calc_bhar_mbh(),
  // except that now we're looking at fixed stellar mass,
  // so we need to calculate the average and median BH
  // mass first.
  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);
  double med_to_avg_bh = 0.5 * ((s2 * s2) * M_LN10);


  if (mbh_med + 3 * s2 < 5) return 0;

  double sm[M_BINS]={0}, t[M_BINS]={0}, mbh_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    mbh_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tmbh_avg = mbh_avg[i] + f*(mbh_avg[i+1]-mbh_avg[i]);
    if (!isfinite(tmbh_avg) || tmbh_avg < 5) continue;
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);

    // Assuming different BHs share the same Eddington ratio
    // distribution at fixed ***halo mass***, the BH accretion
    // rates are enhanced by the same amount that the BH mass
    // is relative to the average BH mass for this halo mass bin.
    tbhar *= exp10fc(mbh_med + med_to_avg_bh - tmbh_avg);
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
  }

  bhar_nd /= nd;
  if (bhar_nd > 1e-15) return bhar_nd;
  return 1e-15;
}


// calculate the average BHAR/SFR ratio as a function of BH mass (Mbh)
// and redshift (z). This function is actually more complicated 
// than the calculation of BHAR, BHMR, and BHERs, because there's a 
// correlation between the SFR and the average stellar mass. 
// Consequently we need to trace the galaxy mass distribution 
// in each halo mass, ***given the BH mass*** that we're looking for.
// Most comments are omitted because they are already written in 
// calc_bhar_mbh().
double calc_bhar_sfr_mbh(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = fabs((1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter));
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1*s1+s2*s2;
  double scatter_tot = sqrt(s3);
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;
  double sfr_nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  double sm[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (scatter * M_LN10 * scatter * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double sfr_fine_bin = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(m - tbhm_avg);

    double tsfr = 0;
    double tnorm = 0;

    // The calculation of SFR at fixed BH mass should 
    // be done carefully, by using the prob calculated in
    // the loop below.
    for (double sm_tmp=3; sm_tmp<=13; sm_tmp+=0.05)
    {
      double bm_tmp = bulge_mass(sm_tmp + mu, 1/(1+z));
      double dbhm_tmp = (m - calc_bh_at_bm(bm_tmp, steps[step].smhm)) / s2;
      double dsm_tmp = (sm_tmp - tsm) / scatter;
      double sfr_tmp = sfr_fine_bin * pow(exp10(sm_tmp) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);
      double prob = scatter_tot / (s2 * scatter) * exp(-0.5*dbhm_tmp*dbhm_tmp - 0.5*dsm_tmp*dsm_tmp);
      tnorm += prob;
      tsfr += sfr_tmp * prob;

    }
    tsfr /= tnorm;


    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    if (isfinite(tbhar)) bhar_nd += weight*tnd*tbhar;
    if (isfinite(tsfr)) sfr_nd += weight*tnd*tsfr;
  }

  double bhar_sfr = bhar_nd / sfr_nd;
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}



// calculate the average BHAR/SFR ratio as a function of Mstar and z.
// This function is actually more complicated than the calculation of
// BHAR, BHMR, and BHERs, because there's a correlation between the 
// SFR and the average stellar mass. Consequently we need to trace
// the galaxy mass distribution in each halo mass, ***given the
// BH mass*** that we're looking for.
// Most comments are omitted because they are already written in 
// calc_bhar_mbh() or calc_bhar_sfr_mbh().
double calc_bhar_sfr_mstar(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = fabs((1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter));
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double s3 = s1 * s1;

  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;
  double sfr_nd = 0;

  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);
  double med_to_avg_bh = 0.5 * ((s2 * s2) * M_LN10);

  if (mbh_med + 3 * s2 < 5) return 0;

  double sm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (s1 * M_LN10 * s1 * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(mbh_med + med_to_avg_bh - tbhm_avg);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);

    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    if (isfinite(tbhar)) bhar_nd += weight*tnd*tbhar;
    if (isfinite(tsfr)) sfr_nd += weight*tnd*tsfr;
  }
  double bhar_sfr = bhar_nd / sfr_nd;
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}

// calculate the average SFR as a function of stellar mass (Mstar)
// and redshift (z). Most comments are omitted because they are 
// already written in calc_bhar_mbh().
double calc_sfr_mstar(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1 * s1;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double sfr_nd = 0;

  double sm[M_BINS]={0}, t[M_BINS]={0};
  double sfr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (s1 * M_LN10 * s1 * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    sfr_nd += weight*tnd*tsfr;
  }
  sfr_nd /= nd;
  
  if (sfr_nd > 1e-15) return sfr_nd;
  return 1e-15;
}




// calculate the average SBHAR/SSFR ratio as a function of black hole mass (Mbh)
// and redshift (z). This function is actually more complicated than the calculation of
// BHAR, BHMR, and BHERs, because there's a correlation between the 
// SFR and the average stellar mass. Consequently we need to trace
// the galaxy mass distribution in each halo mass, ***given the
// BH mass*** that we're looking for.
double calc_sbhar_ssfr_mbh(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = fabs((1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter));
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1*s1+s2*s2;
  double scatter_tot = sqrt(s3);
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;
  double sfr_nd = 0;
  double sm_avg_nd = 0;


  double bhm[M_BINS]={0}, bhm_avg[M_BINS]={0}, t[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  double sm[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (scatter * M_LN10 * scatter * M_LN10));
  double med_to_avg_bh = exp(0.5 * (s3 * M_LN10 * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    if (sm[i+1] == 0) f=0;
    if (sm[i] == 0) continue;
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    tbhar *= exp10fc(m - tbhm_avg);
    double sfr_fine_bin = sfr[i] + f*(sfr[i+1]-sfr[i]);

    double tsm_avg = 0;
    double tsfr = 0;
    double tnorm = 0;
    for (double sm_tmp=0; sm_tmp<=15; sm_tmp+=0.05)
    {
      double bm_tmp = bulge_mass(sm_tmp + mu, 1/(1+z));
      double dbhm_tmp = (m - calc_bh_at_bm(bm_tmp, steps[step].smhm)) / s2;
      double dsm_tmp = (sm_tmp - tsm) / scatter;
      double sfr_tmp = sfr_fine_bin * pow(exp10(sm_tmp) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);
      double prob = scatter_tot / (s2 * scatter) * exp(-0.5*dbhm_tmp*dbhm_tmp - 0.5*dsm_tmp*dsm_tmp);
      tnorm += prob;
      tsfr += sfr_tmp * prob;
      tsm_avg += exp10(sm_tmp) * med_to_avg_star * prob;
    }
    tsm_avg /= tnorm;
    tsfr /= tnorm;
    
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    if (isfinite(tbhar)) bhar_nd += weight*tnd*tbhar;
    if (isfinite(tsfr)) sfr_nd += weight*tnd*tsfr;
    if (isfinite(tsm_avg)) sm_avg_nd += weight * tnd * tsm_avg;
  }
  bhar_nd /= nd; //divide by the number density
  sfr_nd /= nd;
  sm_avg_nd /= nd;
  double bhar_sfr = bhar_nd / sfr_nd / (exp10(m)) * sm_avg_nd;
  double bm_scaling = (m - steps[step].smhm.bh_beta) / steps[step].smhm.bh_gamma + 11;
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}

// calculate the average specific BHAR (BHAR / BH mass) as a 
// function of stellar mass (Mstar) and redshift (z).
double calc_sbhar_mstar(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = fabs((1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter));
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1 * s1;
  if (!s3) 
  {
    // fprintf(stderr, "s3=%f, s1=%f, scatter=%f, gamma=\n");
    return 1e-15;
  }
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;
  double bhm_avg_nd = 0;
  double sfr_nd = 0;
  double sm_avg_nd = 0;

  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);
  double med_to_avg_bh = (0.5 * ((s2 * s2) * M_LN10));

  if (mbh_med + 3*s2 < 5) return 0; 
  double sm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  double bhm[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (scatter * M_LN10 * scatter * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tbm = bulge_mass(m + mu, 1/(1+z));
    double tbhm = calc_bh_at_bm(tbm, steps[step].smhm);
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(mbh_med + med_to_avg_bh - tbhm_avg);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);


    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    if (isfinite(tbhar)) bhar_nd += weight*tnd*tbhar;
    if (isfinite(tsfr)) sfr_nd += weight*tnd*tsfr;
    if (isfinite(tbhm)) bhm_avg_nd += weight*tnd*exp10(tbhm +  med_to_avg_bh);
  }

  bhar_nd /= nd;
  sfr_nd /= nd;
  bhm_avg_nd /= nd;

  double sbhar = bhar_nd / bhm_avg_nd;
  double tbm = bulge_mass(m + mu, 1/(1+z));
  double tbhm = calc_bh_at_bm(tbm, steps[step].smhm);
  if (sbhar > 1e-15) return sbhar;
  return 1e-15;
}




// calculate the average SBHAR/SSFR ratio as a function of Mstar and z.
// This function is actually more complicated than the calculation of
// BHAR, BHMR, and BHERs, because there's a correlation between the 
// SFR and the average stellar mass. Consequently we need to trace
// the galaxy mass distribution in each halo mass, ***given the
// BH mass*** that we're looking for.
double calc_sbhar_ssfr_mstar(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = fabs((1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter));
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1 * s1;
  if (!s3) 
  {
    // fprintf(stderr, "s3=%f, s1=%f, scatter=%f, gamma=\n");
    return 1e-15;
  }
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;
  double bhm_avg_nd = 0;
  double sfr_nd = 0;
  double sm_avg_nd = 0;

  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);
  double med_to_avg_bh = (0.5 * ((s2 * s2) * M_LN10));

  if (mbh_med + 3*s2 < 5) return 0;

  double sm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  double bhm[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (scatter * M_LN10 * scatter * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tbm = bulge_mass(m + mu, 1/(1+z));
    double tbhm = calc_bh_at_bm(tbm, steps[step].smhm);
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(mbh_med + med_to_avg_bh - tbhm_avg);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);


    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    if (isfinite(tbhar)) bhar_nd += weight*tnd*tbhar;
    if (isfinite(tsfr)) sfr_nd += weight*tnd*tsfr;
    if (isfinite(tbhm)) bhm_avg_nd += weight*tnd*exp10(tbhm) * exp10fc(med_to_avg_bh);
  }

  bhar_nd /= nd;
  sfr_nd /= nd;
  bhm_avg_nd /= nd;
  double bhar_sfr = bhar_nd / sfr_nd / (bhm_avg_nd / (exp10(m)));
  double tbm = bulge_mass(m + mu, 1/(1+z));
  double tbhm = calc_bh_at_bm(tbm, steps[step].smhm);
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}



// calculate the average BHER as a function of Mbh and z.
double calc_bher_mbh(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bher_nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0};
  double bher[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bher[i] = 4.5e8 * ((1.0-f)*(steps[step].bh_acc_rate[i] / steps[step].bh_mass_avg[i] * steps[step].smhm.bh_efficiency_rad) + 
                  f*steps[step+1].bh_acc_rate[i] / steps[step+1].bh_mass_avg[i] * steps[step+1].smhm.bh_efficiency_rad);
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbher = bher[i] + f*(bher[i+1]-bher[i]);
    if (!isfinite(tbher)) continue;
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    bher_nd += weight*tnd*tbher;
  }
  bher_nd /= nd;
  
  if (bher_nd > 1e-15) return bher_nd;
  return 1e-15;
}

// calculate the average BHER as a function of Mstar and z.
double calc_bher_mstar(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step].smhm.scatter);
  double s3 = s1*s1;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bher_nd = 0;

  double sm[M_BINS]={0}, t[M_BINS]={0};
  double bher[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    bher[i] = 4.5e8 * ((1.0-f)*(steps[step].bh_acc_rate[i] / steps[step].bh_mass_avg[i] * steps[step].smhm.bh_efficiency_rad) + 
                  f*steps[step+1].bh_acc_rate[i] / steps[step+1].bh_mass_avg[i] * steps[step+1].smhm.bh_efficiency_rad);
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbher = bher[i] + f*(bher[i+1]-bher[i]);
    if (!isfinite(tbher)) continue;
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    bher_nd += weight*tnd*tbher;
  }
  bher_nd /= nd;
  if (bher_nd > 1e-15) return bher_nd;
  return 1e-15;
}



// calculate the average BHMR as a function of Mbh and z.
double calc_bhmr_mbh(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhmr_nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhmr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) 
  {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bhmr[i] = (1.0-f)*steps[step].bh_merge_rate[i] + f*steps[step+1].bh_merge_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhmr = bhmr[i] + f*(bhmr[i+1]-bhmr[i]);
    tbhmr *= exp10fc(m - tbhm_avg);
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    nd += weight*tnd;
    bhmr_nd += weight*tnd*tbhmr;
  }
  bhmr_nd /= nd;
  if (bhmr_nd > 1e-15) return bhmr_nd;
  return 1e-15;
}

// calculate the average BHMR as a function of Mstar and z.
double calc_bhmr_mstar(double m, double z) 
{
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = fabs((1.0-f)*(steps[step].smhm.scatter) + f*(steps[step].smhm.scatter));
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double s3 = s1*s1;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhmr_nd = 0;

  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);
  double med_to_avg_bh = (0.5 * ((s2 * s2) * M_LN10));

  if (mbh_med + 3*s2 < 5) return 0;

  double sm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhmr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    bhmr[i] = (1.0-f)*steps[step].bh_merge_rate[i] + f*steps[step+1].bh_merge_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) 
  {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }

    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg) || tbhm_avg < 5) continue;
    double tbhmr = bhmr[i] + f*(bhmr[i+1]-bhmr[i]);
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    tbhmr *= exp10fc(mbh_med + med_to_avg_bh - tbhm_avg);
    nd += weight*tnd;
    bhmr_nd += weight*tnd*tbhmr;
  }
  bhmr_nd /= nd;
  if (bhmr_nd > 1e-15) return bhmr_nd;
  return 1e-15;
}

// Calculate the cosmic black hole mass density at a given
// redshift (z), and a threshold in Eddington ratio (thresh_ledd)
// or in bolometric luminosity (thresh_lbol). thresh_ledd
// and thresh_lbol should be set to zero if no limit is applied.
double cosmic_bh_density(double z, double thresh_ledd, double thresh_lbol, struct smf_fit *fit) 
{
  int64_t i, step;
  double f;

  // Find out which snapshot to look at
  calc_step_at_z(z, &step, &f);
  double mnd = 0;

  // Alpha and delta will be useful when a threshold on Eddington ratio or
  // bolometric luminosity is applied.
  double alpha = (1.0-f)*steps[step].smhm.bh_alpha + f*(steps[step+1].smhm.bh_alpha);
  double delta = (1.0-f)*steps[step].smhm.bh_delta + f*(steps[step+1].smhm.bh_delta);
  double bh_eta_crit = (1.0-f)*steps[step].smhm.bh_eta_crit + f*(steps[step+1].smhm.bh_eta_crit);
  for (i=0; i<M_BINS; i++) 
  {
    // mass- and redshift-dependent AGN duty cycles.
    double dc = (1.0-f)*steps[step].smhm.bh_duty + f*(steps[step+1].smhm.bh_duty);
    double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    double sn = dc/doublePL_frac_above_thresh(BHER_EFF_MIN, alpha, delta, 1.0, fit);

    // average BH mass for each halo mass bin, interpolated between two snapshots.
    double bhm = (1.0-f)*steps[step].bh_mass_avg[i] + f*(steps[step+1].bh_mass_avg[i]);
    // halo number densities, interpolated between two snapshots.
    double nd = (1.0-f)*steps[step].t[i] + f*(steps[step+1].t[i]);

    // If the threshold is in the form of bolometric luminosity,
    // convert it into the one in Eddington ratio.
    if (thresh_lbol) 
    {
      double lbhm = (1.0-f)*steps[step].log_bh_mass[i] + f*(steps[step+1].log_bh_mass[i]);
      double lum_i = 72.5 - 2.5*(thresh_lbol - 7);
      thresh_ledd = -(lum_i + 5.26)/2.5 - lbhm;

      // If we assume a non-linear relation between the radiative and total
      // Eddington ratios, we also need to convert the threshold (radiative) 
      // Eddington ratio into the total one.
      if (nonlinear_luminosity) 
      {
        const double log_of_2 = M_LN2 / M_LN10;
        if (thresh_ledd > log_of_2) 
        {
          double bher_norm = exp10(thresh_ledd - log_of_2);
          if (!isfinite(bher_norm)) continue;
          thresh_ledd = (bher_norm - 1.0 + M_LN2)/M_LN10;
        }
        else if (thresh_ledd < bh_eta_crit) 
        {
          thresh_ledd = (thresh_ledd + bh_eta_crit)/2.0;
        }
      }
    }

    // If the threshold is in Eddington ratio, then just calculate the fraction of
    // BHs that lie above the threshold.
    if (thresh_ledd) 
    {
      double beta = (1.0-f)*steps[step].bh_eta[i] + f*(steps[step+1].bh_eta[i]);
      double t_ef = thresh_ledd - beta;
      double f_active = doublePL_frac_above_thresh(t_ef, alpha, delta, sn, fit);
      mnd += nd*bhm*f_active;
    } 

    // If no threshold is applied, then add everything up.
    else 
    {
      mnd += nd*bhm;
    }
  }
  return mnd;
}

// Calculate cosmic BH mass density contributed by halos between
// [Mh_low, Mh_high], at a given redshift (z).
double cosmic_bh_density_split(double z, double Mh_low, double Mh_high, struct smf_fit *fit) 
{
  int64_t i, step;
  double f;
  double mh;
  calc_step_at_z(z, &step, &f);
  double mnd = 0;
  double nd_tot = 0;
  
  for (i=0; i<M_BINS; i++) 
  {
    mh = M_MIN + (i + 0.5) * INV_BPDEX;
    if (mh < Mh_low || mh > Mh_high) continue;
    double bhm = exp10((1.0-f)*steps[step].log_bh_mass[i] + f*(steps[step+1].log_bh_mass[i]));
    double nd = (1.0-f)*steps[step].t[i] + f*(steps[step+1].t[i]);
    mnd += nd*bhm;
    nd_tot += nd;
  }
  mnd /= nd_tot;
  return mnd;
}

// Calculate quasar luminosity functions at a given luminosity (Mi)
// and redshift (z).
double calc_quasar_lf_new(double Mi, double z) 
{
  // The input luminosity is in the i band magnitude at z=2,
  // which can be converted in to luminosity (erg/s) like this:
  double lbol = 36 - 0.4 * Mi;
  int64_t step; double sf = 0;
  calc_step_at_z(z, &(step), &sf);
  int64_t i=0;
  double mbh_min = steps[step].bh_mass_min; double mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  if (lbol < LBOL_MIN || lbol > LBOL_MAX) return 0;
  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  // In calc_bh_lum_distribution_full() in calc_sfh.c,
  // we have already calculated the quasar luminosity
  // functions for each BH mass bin. To calculate the
  // total quasar luminosity, we simply add them up.
  double ld = 0;
  for (i=0; i<MBH_BINS; i++)
  {
    double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
    int64_t lbol_b = lbol_f;
    lbol_f -= lbol_b;
    double p1 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b];
    double p2 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b+1];
    ld += p1 + lbol_f * (p2 - p1);
  }
  return norm*ld;
}

// Calculate quasar luminosity functions at a given luminosity (Mi)
// and redshift (z) that are contributed by galaxies with masses
// between [mstar_low, mstar_high].
double calc_quasar_lf_mstar(double Mi, double z, double mstar_low, double mstar_high) 
{
  double lbol = 36 - 0.4 * Mi;
  int64_t step; double sf = 0;
  calc_step_at_z(z, &(step), &sf);
  int64_t i, j, k;
  double mbh_min = steps[step].bh_mass_min; double mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  double bh_scatter = steps[step].smhm.bh_scatter;
  double scatter = steps[step].smhm.scatter;
  double mu = steps[step].smhm.mu;

  int64_t mstar_bpdex = 20;
  double mstar_inv_bpdex = 1.0 / mstar_bpdex;
  if (mstar_high - mstar_low < mstar_inv_bpdex)
  {
    fprintf(stderr, "The (mstar_low, mstar_high) interval is too small. The minimum width is %.4f\n", mstar_inv_bpdex);
    return 1e-19;
  }
  int64_t mstar_bins = (mstar_high - mstar_low) * mstar_bpdex;
  double gauss_norm_bh = 1.0 / sqrt(2*M_PI) / bh_scatter;
  double gauss_norm_gal = 1.0 / sqrt(2*M_PI) / scatter;

  double *lum_func_mstar = malloc(sizeof(double)*mstar_bins);
  memset(lum_func_mstar, 0, sizeof(double)*mstar_bins);

  if (lbol < LBOL_MIN || lbol > LBOL_MAX) return 0;
  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  double ld = 0;
  double qlf_mstar = 0;

  double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
  int64_t lbol_b = lbol_f;
  lbol_f -= lbol_b;
  
  for (i=0; i<mstar_bins; i++)
  {
    double sm = mstar_low + (i + 0.5) * mstar_inv_bpdex;
    double bm = bulge_mass(sm + mu, steps[step].scale);
    double mbh_med = calc_bh_at_bm(bm, steps[step].smhm);
    

    for (j=0; j<M_BINS; j++)
    {
      
      double dsm = (sm - steps[step].log_sm[j]) / scatter;
      double prob_mstar = gauss_norm_gal * exp(-0.5*dsm*dsm) * mstar_inv_bpdex;

      double dc = steps[step].smhm.bh_duty;
      double f_mass = exp((log10(steps[step].bh_mass_avg[j]) - steps[step].smhm.dc_mbh) / 
                        steps[step].smhm.dc_mbh_w);
      f_mass /= (1 + f_mass);
      dc *= f_mass;
      if (dc < 1e-4) dc = 1e-4;


      for (k=0; k<MBH_BINS; k++)
      {
        double mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
        double dmbh = (mbh - mbh_med) / bh_scatter;
        double prob_mbh = gauss_norm_bh * exp(-0.5*dmbh*dmbh) * mbh_inv_bpdex;

        double eta_frac = lbol - 38.1 - mbh - steps[step].bh_eta[j];
        double bher_f = (eta_frac-steps[step].ledd_min[j])*steps[step].ledd_bpdex[j];
        int64_t bher_b = bher_f;
        bher_f -= bher_b;
        if (bher_b < 0) {bher_b = 0; bher_f = 0;}
        else if (bher_b >= BHER_BINS - 1) {bher_b = BHER_BINS - 2; bher_f = 1;}
        double p1 = steps[step].bher_dist[j*BHER_BINS+bher_b];
        double p2 = steps[step].bher_dist[j*BHER_BINS+bher_b+1];
        lum_func_mstar[i] += prob_mstar * prob_mbh * (p1 + bher_f * (p2 - p1)) * steps[step].t[j] * dc;
      }

    }

    qlf_mstar += lum_func_mstar[i];
    // printf("sm=%.3f, smf[%d]/smf_norm=%.3e, lum_func_mstar[%d]=%.3e, qlf_mstar=%.3e\n", sm, 
    //   i, smf[i] / smf_norm, i, lum_func_mstar[i], qlf_mstar);
  }
  return qlf_mstar / 2.5;
}

// Calculate quasar luminosity functions at a given luminosity (Mi)
// and redshift (z) that are contributed by BHs with masses
// between [mbh_low, mbh_high].
double calc_quasar_lf_mbh(double Mi, double z, double mbh_low, double mbh_high) 
{
  double lbol = 36 - 0.4 * Mi;
  int64_t step; double sf = 0;
  calc_step_at_z(z, &(step), &sf);
  int64_t i=0;
  double mbh_min = steps[step].bh_mass_min; double mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  // Find out which Mbh bins we should look at.
  double mbh_low_f = (mbh_low - mbh_min) / mbh_inv_bpdex;
  double mbh_high_f = (mbh_high - mbh_min) / mbh_inv_bpdex;
  int64_t mbh_low_b = mbh_low_f;
  mbh_low_f -= mbh_low_b;
  int64_t mbh_high_b = mbh_high_f;
  mbh_high_f -= mbh_high_b;

  if (lbol < LBOL_MIN || lbol > LBOL_MAX) return 0;
  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  double ld = 0;
  for (i=mbh_low_b; i<=mbh_high_b; i++)
  {
    // This w variable is to account for the incomplete contributions
    // from the first and the last Mbh bins, because generally the mass
    // cuts here is not aligned with the Mbh bins.
    double w = 1; //If the whole bin is included, then w = 1;
    if (i == mbh_low_b) w = 1 - mbh_low_f; //At the smallest Mbh bin, only the
                                          //BHs above the cut will contribute
    else if (i == mbh_high_b) w = mbh_high_f; //The opposite is the case for the 
                                              //largest Mbh bin.
    if (mbh_high_b == mbh_low_b) w = mbh_high_f - mbh_low_f; //If the slice in Mbh
                                                //is too thin to cover a single Mbh bin...


    double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
    int64_t lbol_b = lbol_f;
    lbol_f -= lbol_b;
    double p1 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b];
    double p2 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b+1];
    ld += w * (p1 + lbol_f * (p2 - p1));
  }
  return norm*ld;
}

// To calculate the QLF in slices of Eddington ratios, we simply convert the eta limits
// to Mbh limits. After the conversion, the calculation is the same as in
// calc_quasar_lf_mbh().
double calc_quasar_lf_eta(double Mi, double z, double eta_low, double eta_high) 
{
  double lbol = 36 - 0.4 * Mi;
  int64_t step; double sf = 0;
  calc_step_at_z(z, &(step), &sf);
  int64_t i=0;
  double mbh_min = steps[step].bh_mass_min; double mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  double mbh_low = lbol - 38.1 - eta_high;
  double mbh_high = lbol - 38.1 - eta_low;
  if (mbh_low < mbh_min) mbh_low = mbh_min;
  if (mbh_high > mbh_max) mbh_high = mbh_max;
  // Find out which Mbh bins we should look at.
  double mbh_low_f = (mbh_low - mbh_min) / mbh_inv_bpdex;
  double mbh_high_f = (mbh_high - mbh_min) / mbh_inv_bpdex;
  int64_t mbh_low_b = mbh_low_f;
  mbh_low_f -= mbh_low_b;
  int64_t mbh_high_b = mbh_high_f;
  mbh_high_f -= mbh_high_b;

  if (lbol < LBOL_MIN || lbol > LBOL_MAX) return 0;
  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  double ld = 0;
  for (i=mbh_low_b; i<=mbh_high_b; i++)
  {
    // This w variable is to account for the incomplete contributions
    // from the first and the last Mbh bins, because generally the mass
    // cuts here is not aligned with the Mbh bins.
    double w = 1; //If the whole bin is included, then w = 1;
    if (i == mbh_low_b) w = 1 - mbh_low_f; //At the smallest Mbh bin, only the
                                          //BHs above the cut will contribute
    else if (i == mbh_high_b) w = mbh_high_f; //The opposite is the case for the 
                                              //largest Mbh bin.
    if (mbh_high_b == mbh_low_b) w = mbh_high_f - mbh_low_f; //If the slice in Mbh
                                                //is too thin to cover a single Mbh bin...


    double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
    //fprintf(stderr, "lbol=%f, LBOL_MIN=%f, LBOL_BPDEX=%f\n", lbol, LBOL_MIN, LBOL_BPDEX);
    int64_t lbol_b = lbol_f;
    lbol_f -= lbol_b;
    double p1 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b];
    double p2 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b+1];
    ld += w * (p1 + lbol_f * (p2 - p1));
    //fprintf(stderr, "step=%d, i=%d, Mbh=%f, p1=%e. p2=%e, lbol_f=%f, lbol_b=%d\n", step, i, mbh_min + (i + 0.5) * mbh_inv_bpdex, p1, p2, lbol_f, lbol_b);
  }
  return norm*ld;
}

// The probability distribution of Eddington ratio distributions at
// a certain ***absolute*** ***radiative*** Eddington ratio (ledd),
// a typical Eddington ratio (bher_char), double power-law indices
// of Eddington ratio distributions (alpha, delta), a provided 
// normalization of Eddington ratio distribution (bh_prob_norm),
// and a critical Eddington ratio value where the scaling between
// the radiative and total Eddington ratios start to change (bh_eta_crit).
// See Section 2.7 of Zhang et al. (2021).
double _prob_of_ledd_nonlinear(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit) 
{
  double bher = ledd;
  double bher_norm = 1;
  const double log_of_2 = M_LN2 / M_LN10;

  // When the radiative Eddington ratio is too high,
  // we adopt a log relation between the radiative
  // and total Eddington ratios to account for the
  // trapped photons at high Eddington ratios.
  if (bher > log_of_2) 
  {
    bher_norm = exp10(ledd - log_of_2);
    if (!isfinite(bher_norm)) return 0;
    bher = (bher_norm - 1.0 + M_LN2)/M_LN10;
  }

  // When the Eddington ratio is too low, switch
  // to another scaling to account for the kinetic
  // energy output of AGNs.
  else if (bher < bh_eta_crit) 
  {
    bher_norm = 0.5;
    bher = (ledd + bh_eta_crit)/2.0;
  }
  bher-=bher_char;
  double retval = bher_norm*bh_prob_norm / (exp10(bher*alpha) + exp10(bher*delta));
  if (!isfinite(retval)) 
  {
    return 0;
  }
  return retval;
}

// The probability distribution of Eddington ratio distributions at
// a certain ***absolute*** ***kinetic*** Eddington ratio (ledd),
// a typical Eddington ratio (bher_char), double power-law indices
// of Eddington ratio distributions (alpha, delta), a provided 
// normalization of Eddington ratio distribution (bh_prob_norm),
// and a critical Eddington ratio value where the scaling between
// the radiative and total Eddington ratios start to change (bh_eta_crit).
// See Section 2.7 of Zhang et al. (2021).
double _prob_of_ledd_kinetic(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit) 
{
  double bher = ledd;
  double bher_norm1 = 1, bher_norm2 = 1;
  // The maximum value that the kinetic ER can reach is 1/4 the critical ER, which is 10**(-1.5) here.
  // Above that there would be no viable total ER solution.
  if (bher > bh_eta_crit - LOG10_4) return 0;

  double edd = exp10(ledd);
  double edd_crit = exp10(bh_eta_crit);


  // Solve for the total Eddington ratio
  double sqrt_delta = sqrt(1 - 4 * edd / edd_crit);

  double edd_tot1 = 0.5 * edd_crit * (1 + sqrt_delta);
  double edd_tot2 = 0.5 * edd_crit * (1 - sqrt_delta);

  double ledd_tot1 = log10(edd_tot1);
  double ledd_tot2 = log10(edd_tot2);

  // And the corresponding Jacobian
  bher_norm1 = fabs((edd_crit - 2 * edd_tot1) / (edd_crit - edd_tot1));
  bher_norm2 = fabs((edd_crit - 2 * edd_tot2) / (edd_crit - edd_tot2));

  // Calculate the ***total*** ER distribution 
  ledd_tot1 -= bher_char;
  ledd_tot2 -= bher_char;

  // and map them to the ***kinetic*** ER distribution
  double retval = bh_prob_norm * (bher_norm1 / (exp10(ledd_tot1*alpha) + exp10(ledd_tot1*delta)) + 
                                  bher_norm2 / (exp10(ledd_tot2*alpha) + exp10(ledd_tot2*delta)));
  if (!isfinite(retval)) 
  {
    return 0;
  }
  return retval;
}

// Similar to _prob_of_ledd_nonlinear(), but assuming the total Eddington
// ratio always equals the radiative one, i.e., no kinetic energy output.
double _prob_of_ledd_linear(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit) 
{
  double bher = ledd;
  double bher_norm = 1;
  bher-=bher_char;
  return bher_norm*bh_prob_norm / (exp10(bher*alpha) + exp10(bher*delta));
}

// Calculate the quasar probability distribution functions at given 
// AGN bolometric luminosity (l), stellar mass (m), and redshift (z).
// NOTE: This function is deprecated for now, because it uses Eddington 
// ratio distributions for each ***halo*** mass bin, whereas in the newest
// implementation, we calculate the luminosity distributions for all 
// ***BH mass*** bins.
double calc_qpdf_at_l_m_z(double l, double m, double z) 
{
  int64_t i;
  int64_t step;
  double f;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &step, &f);
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = steps[step].smhm.scatter;

  double lum_i = 72.5 - 2.5*(l - 7);

  // The observational data are not corrected for Compton-thick
  // obscurations, but our model includes them. So we need to 
  // calculate the correction and apply it to model predictions, 
  // before comparing them with observations.
  double Lx = find_Lx_at_Lbol(l);
  double psi_ctn = psi(Lx, z); // psi_ctn is the fraction of Compton-thin
                              // obscured objects relative to non-obscured
                              // objects, which is assumed to be the fraction
                              // of Compton-thick objects. So in the end
                              // we would divide the result with (1 + psi_ctn)
                              // to get non-obscured & Compton-think obscured
                              // only values.
  double acc_rate_norm = -(lum_i+5.26)/2.5;

  // Get the required BH mass and Eddington ratio to calculate number densities.
  double med_bh_m = calc_bh_at_bm(bulge_mass(m, 1.0/(1.0+z)), steps[step].smhm);
  double med_edd_r = acc_rate_norm - med_bh_m;
  double tw = 0, ld=0;

  // Go over each halo mass bin and collect BHs with masses and Eddington ratios
  // obtained above.
  for (i=0; i<M_BINS; i++) 
  {
    double dm = (steps[step].log_sm[i]+steps[step].smhm.mu - m)/s1;
    double w = exp(-0.5*dm*dm)*steps[step].t[i];
    tw += w;
    double eta_frac = med_edd_r - steps[step].bh_eta[i];
    if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) continue;
    double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
    int64_t bher_b = bher_f;
    bher_f -= bher_b;

    double p1 = steps[step].bher_dist_full[i*BHER_BINS+bher_b];
    double p2 = steps[step].bher_dist_full[i*BHER_BINS+bher_b+1];
    if (bher_b >= BHER_BINS-1) p2 = p1;
    double prob = p1 + bher_f*(p2-p1);

    // mass-dependent modulation of AGN duty cycles
    double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    double dc = steps[step].smhm.bh_duty * f_mass;
    if (dc < 1e-4) dc = 1e-4;
    prob *= dc;
    ld += w*prob;
  }
  if (tw>0) return (ld/tw/(1+psi_ctn)); //correct for Compton-thick obscurations.
  return 1e-15;
}

// Calculate the quasar probability distribution functions at given 
// AGN bolometric luminosity (l), stellar mass (m), and redshift (z).
double calc_qpdf_at_l_m_z_new(double lbol, double m, double z) 
{
  int64_t i;
  int64_t step;
  double f;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &step, &f);
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }

  // Get the scatter around the BHBM relation, and the median BH mass
  // for this particular stellar mass and redshift.
  double bh_scatter = steps[step].smhm.bh_scatter;
  double bm = bulge_mass(m, steps[step].scale);
  double mbh_med = calc_bh_at_bm(bm, steps[step].smhm);
  double gauss_norm = 1 / sqrt(2 * M_PI) / bh_scatter;

  // The observational data are not corrected for Compton-thick
  // obscurations, but our model includes them. So we need to 
  // calculate the correction and apply it to model predictions, 
  // before comparing them with observations.
  double lx = find_Lx_at_Lbol(lbol);
  double psi_ctn = psi(lx, z);

  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  double qpdf = 0;

  // Go over each BH mass bin and count how many BHs in each bin
  // contribute to the BHs hosted by the galaxies with mass m.
  for (i=0; i<MBH_BINS; i++) 
  {
    double mbh = steps[step].bh_mass_min + (i + 0.5) * mbh_inv_bpdex;
    double dmbh = (mbh - mbh_med)/bh_scatter;
    
    // prob_mbh is the fraction of BHs hosted by galaxies of mass m
    // that are of mass mbh.
    double prob_mbh = exp(-0.5*dmbh*dmbh)*gauss_norm * mbh_inv_bpdex;

    // prob_lbol is simply calculated from the luminosity distributions.
    double prob_lbol;
    if (lbol < LBOL_MIN || lbol > LBOL_MAX) 
      prob_lbol = 0;
    else
    {
      double lbol_f = (lbol-LBOL_MIN)*LBOL_BPDEX;
      int64_t lbol_b = lbol_f;
      lbol_f -= lbol_b;


      double p1 = steps[step].lum_dist_full[i*LBOL_BINS+lbol_b];
      double p2 = steps[step].lum_dist_full[i*LBOL_BINS+lbol_b+1];
      if (lbol_b >= LBOL_BINS-1) p2 = p1;
      prob_lbol = p1 + lbol_f*(p2-p1);
    }
    qpdf += prob_mbh * prob_lbol;
  }

  return qpdf / (1 + psi_ctn); //correct for Compton-thick obscurations.
}


// double calc_qpdf_at_sBHAR_m_z(double sBHAR, double m, double z) {
//   int64_t i;
//   int64_t step;
//   double f;
//   calc_step_at_z(z, &step, &f);
//   if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
//   double s1 = steps[step].smhm.scatter;
//   double tw = 0, ld=0;
//   for (i=0; i<M_BINS; i++) 
//   {
//     double dm = (steps[step].log_sm[i]+steps[step].smhm.mu - m)/s1;
//     double w = exp(-0.5*dm*dm)*steps[step].t[i];

//     double Lbol = log10(2.6e35) + m + sBHAR;
//     double lum_i = 72.5 - 2.5*(Lbol - 7);
//     double Lx = find_Lx_at_Lbol(Lbol);
//     double acc_rate_norm = -(lum_i+5.26)/2.5;
//     double med_bh_m = calc_bh_at_bm(bulge_mass(m, 1.0/(1.0+z)), steps[step].smhm);

//     double med_edd_r = acc_rate_norm - med_bh_m + steps[step].smhm.eta_mu;
//     double psi_ctn = psi(Lx, z);

//     tw += w;
//     double eta_frac = med_edd_r - steps[step].bh_eta[i];
    
//     double prob;
//     if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) 
//       prob = 0;
//     else
//     {
//       double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
//       int64_t bher_b = bher_f;
//       bher_f -= bher_b;

//       double p1 = steps[step].bher_dist_full[i*BHER_BINS+bher_b];
//       double p2 = steps[step].bher_dist_full[i*BHER_BINS+bher_b+1];
//       if (bher_b >= BHER_BINS-1) p2 = p1;
//       prob = p1 + bher_f*(p2-p1);
//       prob /= (1 + psi_ctn);
//     }
    
//     double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
//     f_mass = f_mass / (1 + f_mass);
//     double dc = steps[step].smhm.bh_duty * f_mass;
//     if (dc < 1e-4) dc = 1e-4;
//     prob *= dc;
//     ld += w*prob;

//   }
//   if (tw>0) return (ld/tw);
//   return 1e-15;
// }

// Calculate the quasar probability distribution functions at given 
// specific BH accretio nrate (sBHAR), stellar mass (m), and redshift (z).
// The observational data are taken from Aird et al. (2018). For the 
// definition of sBHAR, see Section 2.8 of Zhang et al. (2018) (where
// it is called sLx)
double calc_qpdf_at_sBHAR_m_z_new(double sBHAR, double m, double z) 
{
  int64_t i;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }

  // Get the scatter around the BHBM relation, and the median BH mass
  // for this particular stellar mass and redshift.
  double bh_scatter = steps[step].smhm.bh_scatter;
  double bm = bulge_mass(m, steps[step].scale);
  double mbh_med = calc_bh_at_bm(bm, steps[step].smhm);
  double gauss_norm = 1 / sqrt(2 * M_PI) / bh_scatter;

  // Get the luminosity from the sBHAR. We will use this 
  // to find the probability distribution values.
  double lx = log10(2.6e35 / 25) + m + sBHAR;
  // double lum_i = 72.5 - 2.5*(Lbol - 7);

  // Find out the bolometric luminosity corresponding to the sBHAR
  // values
  double lbol = find_Lbol_at_Lx(lx);

  // The observational data are not corrected for Compton-thick
  // obscurations, but our model includes them. So we need to 
  // calculate the correction and apply it to model predictions, 
  // before comparing them with observations.
  double psi_ctn = psi(lx, z);

  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;
  double eta_mu = steps[step].smhm.eta_mu;
  double qpdf = 0;

  // Go over each BH mass bin and count how many BHs in each bin
  // contribute to the BHs hosted by the galaxies with mass m.
  for (i=0; i<MBH_BINS; i++) 
  {
    double mbh = steps[step].bh_mass_min + (i + 0.5) * mbh_inv_bpdex;
    double dmbh = (mbh - mbh_med)/bh_scatter;
    // prob_mbh is the fraction of BHs hosted by galaxies of mass m
    // that are of mass mbh.
    double prob_mbh = exp(-0.5*dmbh*dmbh)*gauss_norm * mbh_inv_bpdex;

    // prob_lbol is simply calculated from the luminosity distributions.
    double prob_lbol;
    if (lbol + eta_mu < LBOL_MIN || lbol + eta_mu > LBOL_MAX) 
      prob_lbol = 0;
    else
    {
      double lbol_f = (lbol + eta_mu -LBOL_MIN)*LBOL_BPDEX;
      int64_t lbol_b = lbol_f;
      lbol_f -= lbol_b;


      double p1 = steps[step].lum_dist_full[i*LBOL_BINS+lbol_b];
      double p2 = steps[step].lum_dist_full[i*LBOL_BINS+lbol_b+1];
      if (lbol_b >= LBOL_BINS-1) p2 = p1;
      prob_lbol = p1 + lbol_f*(p2-p1);
    }
    qpdf += prob_mbh * prob_lbol;

  }
  return qpdf / (1 + psi_ctn); //correct for Compton-thick obscurations.
}

// Calculate the observed cosmic star formation rates at a
// given redshift (z).
double calc_cosmic_sfr(double z) 
{
  int64_t step;
  double f;
  double sfr;

  // Find out which snapshot to look at.
  calc_step_at_z(z, &(step), &(f));
  sfr = steps[step].observed_cosmic_sfr;

  // interpolate with the next snapshot, if there
  // is a next snapshot.
  if (f>-1) sfr += f*(steps[step+1].observed_cosmic_sfr-sfr);
  if (sfr < 1e-10) sfr = 1e-10;
  return sfr;
}

// Calculate the observed cosmic black hole accretion 
// rate at a given redshift (z).
double calc_cosmic_bhar(double z) 
{
  int64_t step;
  double f;
  double bhar;
  calc_step_at_z(z, &(step), &(f));
  bhar = steps[step].observed_cosmic_bhar;
  if (f>-1) bhar += f*(steps[step+1].observed_cosmic_bhar-bhar);
  if (bhar < 1e-10) bhar = 1e-10;
  return bhar;
}

// The helper function to calculate the chi2 for 
// stellar mass functions (SMFs) and specific star formation
// rates (SSFRs), given a certain comoving volume (Vc).
double chi2_err_helper(double Vc, void *extra_data) 
{
  double m = *((double *)extra_data);

  // Find out the redshift corresponding to the volume
  double z = comoving_volume_to_redshift(Vc);

  // evaluate the SMF and SSFR
  double smf = evaluate_from_step(z, m, 1);
  if (smf < 1e-15) smf = 1e-15;
  return smf;
}

// The helper function to calculate the chi2 for 
// galaxy quenched fractions (QFs), 
// given a certain comoving volume (Vc).
double chi2_err_helper_qf(double Vc, void *extra_data) 
{
  double m = *((double *)extra_data);
  double z = comoving_volume_to_redshift(Vc);
  double qf = evaluate_from_step(z, m, 2);
  return qf;
}

// The helper function to calculate the chi2 for 
// galaxy UV luminosity functions (UVLFs), 
// given a certain comoving volume (Vc).
double chi2_err_helper_uv(double Vc, void *extra_data) 
{
  double uv = *((double *)extra_data);
  double z = comoving_volume_to_redshift(Vc);
  double uvlf = evaluate_from_step_uv(z, uv);
  if (uvlf < 1e-15) uvlf = 1e-15;
  return uvlf;
}

// Calculate the model chi2, given a series of model data values (m_smf),
// and another series of real data values (r_smf).
double model_chi2(struct obs_smf *m_smf, struct obs_smf *r_smf) 
{
  int i;
  double chi2 = 0;
  double err, real_err;
  struct real_smf *real_smf = r_smf->real_smf;
  struct real_smf *model_smf = m_smf->real_smf;

  // iterate over each data point
  for (i=0; i<r_smf->real_smf_points; i++) 
  {
    // fit_error is the tolerance in the absolute difference between the model and data.
    // If the model vs. data difference is smaller than that, we ignore its contribution
    // to the final chi2. This is chosen as empirical models typically cannot fit to 
    // real data to better than this tolerance.
    double fit_error = (r_smf->linear_errors) ? fabs(real_smf[i].val*(exp10fc(FIT_TOLERANCE)-1.0)) : fabs(FIT_TOLERANCE);
    err = (model_smf[i].val-real_smf[i].val);
    if (fabs(err) <= fit_error) continue;
    if (err > 0) err -= fit_error;
    else err += fit_error;

    // If the model vs. data difference is bigger than the upper error bar, then
    // we adopt the upper error bar to scale the difference;
    if (err > real_smf[i].err_h) err/=real_smf[i].err_h;

    // If the model vs. data difference is smaller than the opposite of the 
    // lower error bar, then we adopt the lower error bar to scale the difference;
    else if (-err > real_smf[i].err_l) err/=real_smf[i].err_l;

    // otherwise, we do a linear interpolation to get a "medium" error bar.
    // See Appendix F of Zhang et al. (2021).
    else 
    {
      real_err = real_smf[i].err_l +
    	(err+real_smf[i].err_l)/(real_smf[i].err_h+real_smf[i].err_l) *
    	(real_smf[i].err_h - real_smf[i].err_l);
      err /= real_err;
    }
    if (!isfinite(err)) 
    {
      fprintf(stderr, "Infinite error at z=%f, m=%f, type=%d, model val=%f, real val=%f\n", 
	      0.5*(r_smf->z_low + r_smf->z_high),
	      real_smf[i].mass, r_smf->type, model_smf[i].val, real_smf[i].val);
    }
    chi2 += err*err;
  }
  return chi2;
}

// Calculate the model chi2, given a single model data value (m_smf),
// and another single real data value (r_smf). This is very similar to 
// model_chi2(), except that this function deals with individual
// data points.
double model_chi2_point(struct obs_smf_point *m_smf, struct obs_smf_point *r_smf) 
{
  double chi2 = 0;
  double err, real_err;
  double fit_error = (r_smf->linear_errors) ? fabs(r_smf->val*(pow(10,FIT_TOLERANCE)-1.0)) : fabs(FIT_TOLERANCE);
  err = (m_smf->val-r_smf->val);
  if (fabs(err) <= fit_error) 
  {
    return 0;
  }
  if (err > 0) err -= fit_error;
  else err += fit_error;

  if (err > r_smf->err_h) err/=r_smf->err_h;
  else if (-err > r_smf->err_l) err/=r_smf->err_l;
  else 
  {
    real_err = r_smf->err_l +
      (err+r_smf->err_l)/(r_smf->err_h+r_smf->err_l) *
      (r_smf->err_h - r_smf->err_l);
    err /= real_err;
  }
  if (!isfinite(err)) 
  {
    fprintf(stderr, "Infinite error at z=%f, m=%f, extra=%f, type=%d, model val=%f, real val=%f\n", 
	    0.5*(r_smf->z_low + r_smf->z_high),
	    r_smf->mass, r_smf->extra, r_smf->type, m_smf->val, r_smf->val);
  }
  chi2 += err*err;
  return chi2;
}

// Calculate chi2's for a series of observed data points (cur_smf).
double calc_single_chi2_err(struct obs_smf *cur_smf) 
{
  int i;
  double chi2 = 0, smf_val, epsilon;
  struct real_smf *real_smf = cur_smf->real_smf;
  struct obs_smf model_smf;
  int real_smf_points = cur_smf->real_smf_points;
  double v_high = comoving_volume(cur_smf->z_high);
  double v_low = comoving_volume(cur_smf->z_low);
  double weight = v_high - v_low;
  double m;

  model_smf.real_smf_points = real_smf_points;

  for (i=0; i<real_smf_points; i++) 
  {
    m = real_smf[i].mass;


    // Calculate the model predictions according to the data types.
    if (cur_smf->type == SMF_TYPE) 
    {
      if (cur_smf->z_low != cur_smf->z_high && !no_z_scaling) 
      {
	      if (cur_smf->z_low < 0.2) 
        {
	        smf_val = adaptiveGauss(chi2_err_helper, &m, v_low, v_high,
	 			  PHI_INTEGRAL_PRECISION*0.01, omp_get_thread_num()+20);
	      }
	      else 
        {
      	  epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
      	  smf_val = adaptiveSimpsons(chi2_err_helper, &m,
      				     v_low, v_high, epsilon, 10);
	      }
	      smf_val /= weight;
      }
      else 
      {
	      smf_val = chi2_err_helper(0.5*(v_low+v_high), &m);
      }
    }

    else if (cur_smf->type == QF_TYPE) 
    {
      smf_val = chi2_err_helper_qf(0.5*(v_low+v_high), &m);
    }


    else if (cur_smf->type == SSFR_TYPE) 
    {
      smf_val = calc_ssfr(m, cur_smf->z_low);
    }

    else if (cur_smf->type == UVLF_TYPE)
    {
      smf_val = chi2_err_helper_uv(0.5 * (v_low + v_high), &m);
    }


    else if (cur_smf->type == CSFR_TYPE) 
    {
      smf_val = calc_cosmic_sfr(real_smf[i].mass); //Actually redshift for CSFR, not mass.
    } 

    else if (cur_smf->type == BHMF_TYPEI_TYPE)
    {
      smf_val = calc_bhmf_typeI(m, cur_smf->z_low);
    }

    model_smf.real_smf[i].mass = m;
    model_smf.real_smf[i].val = (cur_smf->linear_errors) ? smf_val : log10(smf_val);
  }
  chi2 = model_chi2(&model_smf, cur_smf);
  return chi2;
}

// Calculate chi2's for a single observed data points (cur_smf).
// This is very similar to calc_single_chi2_err().
double calc_single_chi2_err_point(struct obs_smf_point *cur_smf) 
{
  double smf_val, epsilon;
  struct obs_smf_point model_smf;
  double v_high = comoving_volume(cur_smf->z_high);
  double v_low = comoving_volume(cur_smf->z_low);
  double weight = v_high - v_low;
  double m;
  double median_z = cbrt(0.5*(cur_smf->z_high*cur_smf->z_high*cur_smf->z_high + cur_smf->z_low*cur_smf->z_low*cur_smf->z_low));


  m = cur_smf->mass;
  if (cur_smf->type == SMF_TYPE) 
  {
    if (cur_smf->z_low != cur_smf->z_high) 
    {
    	epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
    	smf_val = adaptiveSimpsons(chi2_err_helper, &m,
    				   v_low, v_high, epsilon, 10);
      smf_val /= weight;
    }
    else 
    {
      smf_val = chi2_err_helper(v_low, &m);
    }
  }
  else if (cur_smf->type == SSFR_TYPE) 
  {
    smf_val = calc_ssfr(m, cur_smf->z_low);
  }
  
  else if (cur_smf->type == CSFR_TYPE) {
    smf_val = calc_cosmic_sfr(cur_smf->mass); //Actually redshift for CSFR, not mass
  } 
  else if (cur_smf->type == QLF_TYPE) 
  {
    smf_val = calc_quasar_lf_new(m, cur_smf->z_low);
  }
  else if (cur_smf->type == QPDF_TYPE) 
  {
    smf_val = calc_qpdf_at_l_m_z_new(cur_smf->extra, cur_smf->mass, median_z);
  }
  else if (cur_smf->type == QPDF_ETA_TYPE)
  {
    smf_val = calc_qpdf_at_sBHAR_m_z_new(cur_smf->extra, cur_smf->mass, median_z);
  }
  else if (cur_smf->type == ABHMF_TYPE) 
  {
    smf_val = calc_active_bhmf(cur_smf->mass, median_z);
  }

  else if (cur_smf->type == QF_TYPE) 
  {
    smf_val = chi2_err_helper_qf(0.5*(v_low + v_high), &m);
  }

  else if (cur_smf->type == UVLF_TYPE)
  {
    smf_val = chi2_err_helper_uv(0.5 * (v_low + v_high), &m);
  }

  else if (cur_smf->type == BHMF_TYPEI_TYPE)
  {
    smf_val = calc_bhmf_typeI(m, cur_smf->z_low);
  }

  model_smf.mass = m;
  model_smf.val = (cur_smf->linear_errors) ? smf_val : log10(smf_val);
  return model_chi2_point(&model_smf, cur_smf);
}

// The helper function to calculate the halo mass function based 
// on the cached simulation data, given a halo peak mass (mass).
double _mf_cache(struct mf_cache *cache, double mass) 
{
  int i;
  double m, r;
  if (mass >= cache->mass_max) 
  {
    return -1000;
  }
  else 
  {
    m = (mass - cache->mass_min)*cache->inv_mass_spacing;
    if (m < 0) r=((1.0-cache->alpha)*(mass-cache->mass_min) + cache->mf[0]);
    else 
    {
      i = m;
      r = (cache->mf[i] + (cache->mf[i+1]-cache->mf[i])*(m-i));
    }
  }
  return r;
}

// Calculate the halo mass function based on the cached simulation data,
// given a certain scale factor (scale = 1 / (1 + z)) and halo peak mass
// (mass).
double mf_cache(double scale, double mass) 
{
  int i;
  double c1, c2;
  if (scale >= all_mf_caches->scale_max)
    return _mf_cache(all_mf_caches->caches+all_mf_caches->scales-1, mass);
  if (scale < all_mf_caches->scale_min)
    return _mf_cache(all_mf_caches->caches, mass);

  //Find scale
  i = (scale - all_mf_caches->scale_min) * all_mf_caches->avg_inv_scale_spacing;
  while (all_mf_caches->caches[i+1].scale < scale) i++;
  while (all_mf_caches->caches[i].scale > scale) i--;
  c1 = _mf_cache(all_mf_caches->caches + i, mass);
  c2 = _mf_cache(all_mf_caches->caches + i + 1, mass);
  return (c1 + (c2 - c1)*(scale - all_mf_caches->caches[i].scale)*
	  all_mf_caches->caches[i].inv_scale_spacing);
}

// The function to load real observational data points.
void load_real_smf(struct obs_smf_point **cur_smf, int64_t *num_points, char *filename)
{
  FILE *smf_data;
  char buffer[1024];

  // different data type identifiers
  char *data_types[] = {"smf", "ssfr", "cosmic sfr", "quasar lf", "quasar luminosity pdf", "active bhmf", "quasar eddington ratio pdf", "Compton-thin quasar lf", "qf", "uvlf", "type I bhmf"};
  
  // error types indicating if the data points and error bars are in
  // linear or log units.
  char *error_types[] = {"log", "linear"};
  struct obs_smf_point next_smf, conv_smf;
  int64_t n, loaded=0;

  // open the data file
  smf_data = fopen(filename, "r");
  if (!smf_data) 
  {
    printf("Couldn't open file %s for reading!\n", filename);
    exit(1);
  }

  next_smf.type = SMF_TYPE;
  next_smf.linear_errors = 0;
  next_smf.z_low = next_smf.z_high = 0;
  while (fgets(buffer, 1023, smf_data)) 
  {
    n = sscanf(buffer, "%lf %lf %lf %lf", &(next_smf.mass),
	       &(next_smf.val), &(next_smf.err_h), &(next_smf.err_l));
    if (buffer[0] == '#') 
    {
      // Get the lower/upper limits in redshifts.
      if (!strncmp(buffer, "#zlow: ", 7)) next_smf.z_low = atof(&(buffer[7]));
      if (!strncmp(buffer, "#zhigh: ", 8)) next_smf.z_high = atof(&(buffer[8]));
      // Get the error type.
      if (!strncmp(buffer, "#errors: ", 9)) 
      {
	      if (!strncasecmp(&(buffer[9]), "Linear", 6)) next_smf.linear_errors = 1;
	      else  next_smf.linear_errors = 0;
      }
      // Get the data type.
      if (!strncmp(buffer, "#type: ", 7)) 
      {
      	if (!strncasecmp(&(buffer[7]), "ssfr", 4)) next_smf.type = SSFR_TYPE;
      	else if (!strncasecmp(&(buffer[7]), "csfr", 4)) next_smf.type = CSFR_TYPE;
        else if (!strncasecmp(&(buffer[7]), "qlf_ctn", 7)) next_smf.type = QLF_CTN_TYPE;
        else if (!strncasecmp(&(buffer[7]), "qlf", 3)) next_smf.type = QLF_TYPE;
        else if (!strncasecmp(&(buffer[7]), "qpdf_eta", 8)) next_smf.type = QPDF_ETA_TYPE;
        else if (!strncasecmp(&(buffer[7]), "qpdf", 4)) next_smf.type = QPDF_TYPE;
        else if (!strncasecmp(&(buffer[7]), "abhmf", 5)) next_smf.type = ABHMF_TYPE;
        else if (!strncasecmp(&(buffer[7]), "qlf_ctn", 7)) next_smf.type = QLF_CTN_TYPE;
        else if (!strncasecmp(&(buffer[7]), "qf", 2)) next_smf.type = QF_TYPE;
        else if (!strncasecmp(&(buffer[7]), "uvlf", 4)) next_smf.type = UVLF_TYPE;
        else if (!strncasecmp(&(buffer[7]), "bhmf_typei", 10)) next_smf.type = BHMF_TYPEI_TYPE;
      	else  next_smf.type = SMF_TYPE;
      }
      continue;
    }
    if (n!=4 || !next_smf.err_l || !next_smf.err_h) 
    {
      if (strlen(buffer) > 1)
	    fprintf(stderr, "Invalid line (skipping): %s\n", buffer);
      continue;
    }

    // For QPDF data, we need an additional column (next_smf.extra) to contain 
    // either luminosity or specific BH accretion rate (sBHAR or sLx, see Section
    // 2.8 of Zhang et al. 2021).
    if (next_smf.type == QPDF_TYPE || next_smf.type == QPDF_ETA_TYPE) 
    {
      n = sscanf(buffer, "%lf %lf %lf %lf %lf", &(next_smf.mass), &(next_smf.extra),
         &(next_smf.val), &(next_smf.err_h), &(next_smf.err_l));

      if (n!=5 || !next_smf.err_l || !next_smf.err_h) 
      {
        if (strlen(buffer) > 1)
          fprintf(stderr, "Invalid line (skipping): %s\n", buffer);
        continue;
      }
    }
    
    //Add in systematic errors from calculations.
    double syst_err = CALCULATION_TOLERANCE;
    conv_smf = next_smf;
    if (conv_smf.linear_errors) 
    {
      if (next_smf.type == QF_TYPE)
      {
        ;
      }
      else
      {
        conv_smf.err_h = log10((conv_smf.err_h + conv_smf.val)/conv_smf.val);
        if (conv_smf.err_l >= conv_smf.val) conv_smf.err_l = 1.0; //1 dex!
        else 
        {
          conv_smf.err_l = -log10((conv_smf.val-conv_smf.err_l)/conv_smf.val);
          if (conv_smf.err_l < 0) conv_smf.err_l = -conv_smf.err_l;
          if (conv_smf.err_l > 1.0) conv_smf.err_l = 1.0;
        }
        conv_smf.val = log10(conv_smf.val);
        conv_smf.linear_errors = 0;
      }
      
    }
    conv_smf.err_h = sqrt(conv_smf.err_h*conv_smf.err_h + syst_err*syst_err);
    conv_smf.err_l = sqrt(conv_smf.err_l*conv_smf.err_l + syst_err*syst_err);
    if (!((*num_points)%100)) 
    {
      *cur_smf = check_realloc(*cur_smf, sizeof(struct obs_smf_point)*
			       ((*num_points)+100), "Allocating SMFs.");
    }

    cur_smf[0][*num_points] = conv_smf;
    *num_points = (*num_points) + 1;
    loaded++;
  }
  fclose(smf_data);
  if (loaded) 
  {
    if (next_smf.z_high > z_max) z_max = next_smf.z_high;
    if (next_smf.z_low < z_min) z_min = next_smf.z_low;
  }

  fprintf(stderr, "#Loaded %"PRId64" points from %s (type: %s; errors: %s)\n", loaded, filename, data_types[next_smf.type], error_types[next_smf.linear_errors]); 
}

// Set up the cache for standard Gaussian distribution.
void setup_psf(int scatter) 
{
  init_gauss_cache();
  use_obs_psf = scatter;
}

// Read a single float number from the file input.
float readfloat(FILE *input) 
{
  char buffer[50] = {0};
  float result;
  fgets(buffer, 50, input);
  sscanf(buffer, "%f", &result);
  return result;
}

// Load the cached halo mass functions from
// a certain file (filename).
void load_mf_cache(char *filename) 
{
  FILE *input;
  int i;
  float omega_m, omega_l, h0;
  struct mf_cache *mf;
  if (!(input = fopen(filename, "r"))) 
  {
    printf("Couldn't open MF cache %s!\n", filename);
    exit(6);
  }

  all_mf_caches = (struct z_mf_cache *)malloc(sizeof(struct z_mf_cache));

  // Read cosmology parameters, initialize the cosmology, and the time table.
  omega_m = all_mf_caches->omega_m = readfloat(input);
  omega_l = all_mf_caches->omega_l = readfloat(input);
  h0 = all_mf_caches->h0 = readfloat(input);
  init_cosmology(omega_m, omega_l, h0);
  init_time_table(omega_m, h0);

  // Read the minimum and maximum scale factors
  all_mf_caches->scale_min = readfloat(input);
  all_mf_caches->scale_max = readfloat(input);
  // Read the number of scale factors (i.e., # of snapshots)
  all_mf_caches->scales = readfloat(input);
  all_mf_caches->avg_inv_scale_spacing = (float)(all_mf_caches->scales-1) / 
    (all_mf_caches->scale_max - all_mf_caches->scale_min);

  if (! (all_mf_caches->scale_min > 0 && all_mf_caches->scale_min < 1.05 &&
	 all_mf_caches->scale_max > all_mf_caches->scale_min &&
	 all_mf_caches->scale_max < 1.05)) 
  {
    printf("Invalid scales in MF cache (%f - %f), expecting 0<a<1.\n",
	   all_mf_caches->scale_min, all_mf_caches->scale_max);
    exit(7);
  }

  // Allocate the space for the cache from each snapshot
  all_mf_caches->caches = (struct mf_cache *)
    malloc(sizeof(struct mf_cache)*all_mf_caches->scales);

  // Read information for every snapshot
  for (i=0; i<all_mf_caches->scales; i++) 
  {
    mf = &(all_mf_caches->caches[i]);
    // minimum/maximum halo masses
    mf->mass_min = readfloat(input);
    mf->mass_max = readfloat(input);
    // number of halo masses
    mf->masses = readfloat(input);
    // scale factor
    mf->scale = readfloat(input);
    // slope of the halo mass function (before the exponential suppression at the massive end)
    mf->alpha = fabs(readfloat(input));
    if (! (mf->scale >= all_mf_caches->scale_min &&
	   mf->scale <= all_mf_caches->scale_max)) 
    {
      printf("Invalid scale in MF cache index %d: %f\n", i, mf->scale);
      exit(7);
    }
    if (! (mf->alpha > 1 && mf->alpha < 3)) 
    {
      printf("Invalid |alpha| in MF cache: %f (Expected 1 - 3)\n", mf->alpha);
      exit(7);
    }

    mf->mf = (float *)malloc(sizeof(float)*mf->masses);
    // Read in the actual mass functions
    fread(mf->mf, sizeof(float), mf->masses, input);
    mf->inv_mass_spacing = (float)(mf->masses-1)/(mf->mass_max - mf->mass_min);
    // The inverse of the characteristic mass of the Schechter halo mass function.
    mf->inv_m0 = log10(((mf->mass_max - mf->mass_min)*(1.0-mf->alpha)
			+ mf->mf[0] - mf->mf[mf->masses-1])/exp10(mf->mass_max));
  }
  for (i=0; i<all_mf_caches->scales-1; i++)
    all_mf_caches->caches[i].inv_scale_spacing = 
      1.0 / (all_mf_caches->caches[i+1].scale - all_mf_caches->caches[i].scale);
  fclose(input);

  gen_exp10cache(); //Required for using MF cache
}
