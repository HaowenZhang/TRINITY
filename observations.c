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

static double doexp10(double x) {
  double a = exp(M_LN10*x);
  return a;
}

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

double find_Lx_at_Lbol(double Lbol)
{
  double BC = 10.83 * pow(exp10(Lbol) / (1e10 * 3.839e33), 0.28) +
                            6.08 * pow(exp10(Lbol) / (1e10 * 3.839e33), -0.02);
  return Lbol - log10(BC);
}

// ref: Schulze+2015
double F_obs(double Lx)
{
  if (Lx < 43) return 0.985;
  return 0.56 + 1.0 / M_PI * atan((43.89 - Lx) / 0.46);
}

int use_obs_psf = 1;
struct z_mf_cache *all_mf_caches = NULL;
double gauss_cache[GAUSS_CACHE_SIZE];

int output_smf = 0;

void init_gauss_cache(void) {
  int i;
  double m;
  for (i=0; i<GAUSS_CACHE_SIZE; i++) {
    m = (double)i*GAUSS_CACHE_STEP + GAUSS_CACHE_MIN;
    gauss_cache[i] = exp(-0.5*m*m)/sqrt(2.0*M_PI);
  }
}

double schechter_norm(double alpha, double lower, double upper, struct smf_fit *fit) {
  lower = exp10fc(lower);
  upper = exp10fc(upper);
  gsl_sf_result rl, rh;
  if (gsl_sf_gamma_inc_e(alpha, lower, &rl) || gsl_sf_gamma_inc_e(alpha, upper, &rh)) {
    //error happened...
    if (fit) {
      INVALIDATE(fit, "Invalid arguments to incomplete gamma function.");
    }
    return 1;
  }
  return((rl.val-rh.val)/log(10.0));
}

double hyperg_gt1(double a, double b, double c, double z, struct smf_fit* fit)
{
  double g_a, g_b, g_c, g_ab, g_ba, g_ca, g_cb;
  double coef1, coef2;
  double F1_1, F1_2;
  gsl_sf_result r[9];
  if (gsl_sf_gamma_e(a, &(r[0])) || gsl_sf_gamma_e(b, &(r[1])) || gsl_sf_gamma_e(c, &(r[2])) ||
        gsl_sf_gamma_e(a - b, &(r[3])) || gsl_sf_gamma_e(b - a, &(r[4])) || gsl_sf_gamma_e(c - a, &(r[5])) ||
        gsl_sf_gamma_e(c - b, &(r[6])))
  {
    if (fit) 
      {
        INVALIDATE(fit, "Invalid arguments to gamma function.");
      }
  }
  g_a = r[0].val; g_b = r[1].val; g_c = r[2].val;
  g_ab = r[3].val; g_ba = r[4].val; g_ca = r[5].val; g_cb = r[6].val;

  if (gsl_sf_hyperg_2F1_e(a, c - b, a - b + 1, 1 / (1 - z), &(r[7])) || 
       gsl_sf_hyperg_2F1_e(b, c - a, b - a + 1, 1 / (1 - z), &(r[8])))
  {
    if (fit) 
      {
        INVALIDATE(fit, "Invalid arguments to Hypergeometric function.");
      }
  }
  coef1 = pow(1 - z, -a) * g_c * g_ba / (g_b * g_ca);
  coef2 = pow(1 - z, -b) * g_c * g_ab / (g_a * g_cb);
  F1_1 = r[7].val; F1_2 = r[8].val;
  return coef1 * F1_1 + coef2 * F1_2;
}

double doublePL_norm(double a, double b, double lower, double upper, struct smf_fit *fit) { 
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


  int n_pts = 200; // # of points in the numerical integral
  double scale = 1.0 / n_pts * (upper - lower);
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
  double scale_a = pow(scale, a);
  double scale_b = pow(scale, b);

  for (int i = 0; i < n_pts; i++)
  {
    
    // printf("x_left=%.3e, x_right=%.3e, x_center=%.3e, scale=%.3e\n", x_left, x_right, x_center, scale);
    val += 1 / (x_center_a + x_center_b) * dx;
    // x_left = x_right;
    // x_right *= scale;
    // x_center *= scale;
    dx *= scale;
    x_center_a *= scale_a;
    x_center_b *= scale_b;
  }
  
  // Also we need a 1/(ln10) in the normalization.
  return val / log(10.0);
  //return val;
  // return((rl.val-rh.val)/log(10.0));
}



double schechter_inv_avg(double alpha, double lower, double upper, struct smf_fit *fit) {
  double norm = schechter_norm(alpha, lower, upper, fit);
  double norm2 = schechter_norm(alpha+1, lower, upper, fit);
  return norm/norm2;
}

double schechter_frac_above_thresh(double thresh, double alpha, double norm, struct smf_fit *fit) {
  if (thresh < BHER_EFF_MIN) thresh = BHER_EFF_MIN;
  //if (nonlinear_luminosity && thresh < (BHER_MIN-2.0)/2.0)
  //thresh = (BHER_MIN-2.0)/2.0;
  if (thresh >= BHER_EFF_MAX) return 0;
  return (schechter_norm(alpha, thresh, BHER_EFF_MAX, fit)*norm);
}

double doublePL_frac_above_thresh(double thresh, double alpha, double delta, double norm, struct smf_fit *fit) {
  if (thresh < BHER_EFF_MIN) thresh = BHER_EFF_MIN;
  //if (nonlinear_luminosity && thresh < (BHER_MIN-2.0)/2.0)
  //thresh = (BHER_MIN-2.0)/2.0;
  if (thresh >= BHER_EFF_MAX) return 0;
  return (doublePL_norm(alpha, delta, thresh, BHER_EFF_MAX, fit)*norm);
}


//inline 
double syst_cache(double dm) {
  return (exp(-0.5*dm*dm)*(0.5*M_2_SQRTPI*M_SQRT1_2));
}


double gen_inv_sigma(double obs_scatter, double scatter) {
  if (!scatter && (!use_obs_psf || !obs_scatter)) return 0;
  double gauss_sigma = sqrt(obs_scatter*obs_scatter + scatter*scatter);
  return 1.0/gauss_sigma;
}

//inline 
double evaluate_psf(double delta_m, double gauss_inv_sigma) {
  if (!gauss_inv_sigma) return (delta_m ? 0 : 1);
  return(gauss_inv_sigma*syst_cache(delta_m*gauss_inv_sigma));
}


//inline 
double _interp_from_sm(double sm, double *list, int64_t n, char *ok) {
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

double _interp_from_sfr(double sm, double *list, int64_t n) {
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

double _interp_from_sfrac(double sm, double *list, int64_t n) {
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

double _interp_from_uv(double uv, double *list, int64_t n, char *ok) {
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

double _interp_from_std_uv(double uv, double *list, int64_t n) {
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
  


double bulge_mass(double sm, double a) {
  double z_mul = 1.0 - 0.5*(1.0-a);
  return sm+log10(z_mul/(1.0+exp(-1.12882*(sm-10.1993))));
}

// double step_integral_helper(double sm, void *extra_data) {
//   struct step_integral_helper_data *ih = 
//     (struct step_integral_helper_data *)extra_data;
//   double delta_m = sm - ih->sm;
//   double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);
//   //New Kappa:
//   double passive_frac = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass))+1);
//   double passive_frac2 = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass-ih->kappa))+1);

//   // double sfrac1 = _interp_from_sm(sm, ih->sfrac_sm, ih->n, ih->ok);
//   // double sfrac2 = _interp_from_sm(sm, ih->sfrac_sm2, ih->n2, ih->ok2);
//   // double psf1 = (1 - sfrac1) * evaluate_psf(delta_m, ih->gauss_inv_sigma) +
//                 // sfrac1 * evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);
//   // double psf2 = (1 - sfrac2) * evaluate_psf(delta_m, ih->gauss_inv_sigma) +
//   //               sfrac1 * evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);


//   psf = psf*passive_frac + 
//     (1.0-passive_frac2)*evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);

//   if (ih->corr)
//     psf *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
//     // psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
//   double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
//   if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
//   if (!isfinite(val)) 
//     {
//       // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
//       val = 1e-17;
//     }
//   // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
//   // return (psf1*val);
//   return (psf*val);
// }


double step_integral_helper(double sm, void *extra_data) {
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;
  double psf1 = evaluate_psf(delta_m, ih->gauss_inv_sigma);
  //New Kappa:
  // double passive_frac = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass))+1);
  // double passive_frac2 = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass-ih->kappa))+1);

  // In the UM parametrization we no longer have different offsets between SF and Q
  // galaxies, so only one PSF need to be evaluated. 
  // double sfrac1 = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  // double psf1 = (1 - sfrac1) * evaluate_psf(delta_m, ih->gauss_inv_sigma) +
  //               sfrac1 * evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);
  // double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);


  // psf = psf*passive_frac + 
  //   (1.0-passive_frac2)*evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);

  if (ih->corr)
    // psf *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
  if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
  if (!isfinite(val)) 
    {
      // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
      val = 1e-17;
    }
  // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
  return (psf1*val);
  // return (psf*val);
}

// The new step_integral_helper for ssfr. This is because we have to separate the
// lists to be interpolated (i.e., SFR*SMF) for SF and Q galaxies.
double step_integral_helper_ssfr(double sm, void *extra_data) {
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;
  // double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);
  //New Kappa:
  // double passive_frac = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass))+1);
  // double passive_frac2 = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass-ih->kappa))+1);

  double sfrac = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  // double sfrac2 = _interp_from_sm(sm, ih->sfrac_sm2, ih->n2, ih->ok2);
  double psf1 = (1 - sfrac) * evaluate_psf(delta_m, ih->gauss_inv_sigma);
  double psf2 = psf1 / (1 - sfrac) * sfrac;
  // double psf2 = sfrac * evaluate_psf(delta_m, ih->gauss_inv_sigma);


  // psf = psf*passive_frac + 
  //   (1.0-passive_frac2)*evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);

  if (ih->corr)
  {
    // psf *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf2 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  }
  double val1 = _interp_from_sfr(sm, ih->sfr_sm_q, ih->n);
  double val2 = _interp_from_sfr(sm, ih->sfr_sm_sf, ih->n);
  if (ih->f > -1) val1 += ih->f*(_interp_from_sfr(sm, ih->sfr_sm_q2, ih->n2)-val1);
  if (ih->f > -1) val2 += ih->f*(_interp_from_sfr(sm, ih->sfr_sm_sf2, ih->n2)-val2);
  if ((!isfinite(val1)) || (!isfinite(val2))) 
    {
      // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
      val1 = val2 = 1e-17;
    }
  // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
  return (psf1*val1 + psf2*val2);
  // return (psf*val);
}

double step_integral_helper_uv(double uv, void *extra_data) {
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
      // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
      val = 1e-17;
    }
  // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
  //fprintf(stderr, "n_snap=%d, uv=%f, target uv=%f, std_uv=%f, psf=%e, val=%e, val1=%e, val2=%e, f=%f, psf*val=%e\n", ih->n, uv, ih->sm, std_uv, psf, val, val1, val2, ih->f, psf*val);
  return (psf*val);
  // return (psf*val);
}

double step_integral_helper_Q(double sm, void *extra_data) {
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;
  double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);
  //New Kappa:
  // double passive_frac = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass))+1);
  // double passive_frac2 = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass-ih->kappa))+1);

  double sfrac1 = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  // double sfrac2 = _interp_from_sm(sm, ih->sfrac_sm2, ih->n2, ih->ok2);
  double psf1 = (1 - sfrac1) * evaluate_psf(delta_m, ih->gauss_inv_sigma);
  // double psf2 = (1 - sfrac2) * evaluate_psf(delta_m, ih->gauss_inv_sigma) +
  //               sfrac1 * evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);


  // psf = psf*passive_frac + 
  //   (1.0-passive_frac2)*evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);

  if (ih->corr)
    // psf *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
  if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
  if (!isfinite(val)) 
    {
      // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
      val = 1e-17;
    }
  // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
  // fprintf(stderr, "sm=%f, sfrac=%e, psf=%f, val=%f\n", sm, sfrac1, psf1, val);
  return (psf1*val);
  // return (psf*val);
}


double step_integral_helper_SF(double sm, void *extra_data) {
  struct step_integral_helper_data *ih = 
    (struct step_integral_helper_data *)extra_data;
  double delta_m = sm - ih->sm;
  double psf = evaluate_psf(delta_m, ih->gauss_inv_sigma);
  //New Kappa:
  // double passive_frac = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass))+1);
  // double passive_frac2 = 1.0/(exp10fc(-1.3*(sm-ih->passive_mass-ih->kappa))+1);

  double sfrac1 = _interp_from_sfrac(sm, ih->sfrac_sm, ih->n);
  // double sfrac2 = _interp_from_sm(sm, ih->sfrac_sm2, ih->n2, ih->ok2);
  double psf1 = sfrac1 * evaluate_psf(delta_m, ih->gauss_inv_sigma);
  // double psf1 = sfrac1 * evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);
  // double psf2 = (1 - sfrac2) * evaluate_psf(delta_m, ih->gauss_inv_sigma) +
  //               sfrac1 * evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);


  // psf = psf*passive_frac + 
  //   (1.0-passive_frac2)*evaluate_psf(delta_m-ih->kappa, ih->gauss_inv_sigma);

  if (ih->corr)
    // psf *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
    psf1 *= exp10fc(-ih->corr*(delta_m-ih->s_corr));
  double val = _interp_from_sm(sm, ih->list, ih->n, ih->ok);
  if (ih->f > -1) val += ih->f*(_interp_from_sm(sm, ih->list2, ih->n2, ih->ok2)-val);
  if (!isfinite(val)) 
    {
      // printf("sm=%f, delta_m=%f, psf=%f, val=%f\n", sm, delta_m, psf, val);
      val = 1e-17;
    }
  // fprintf(stderr, "sm=%f, sfrac=%e, psf=%f, val=%f\n", sm, sfrac1, psf1, val);
  return (psf1*val);
  // return (psf*val);
}


double evaluate_from_step(double z, double sm, int smf) {
  double phi_corr=1, sm_true, sm_max;
  struct step_integral_helper_data ih;
  calc_step_at_z(z, &(ih.step), &(ih.f));
  struct smf smhm = steps[ih.step].smhm;
  double scatter_corr = smhm.scatter*smhm.scatter*(0.5*M_LN10);

  sm_true = (sm - smhm.mu); // + smhm.kappa*SLOPE_SYSTEMATIC_MASS)/(1.0+smhm.kappa);
  if (!(sm_true >= 0 && sm_true < 20)) {
    fprintf(stderr, "Gah! %f %f %f %f %f %f %f %f %f %f %f\n", sm_true, sm, smhm.sm_0, smhm.v_1, smhm.delta, smhm.beta, smhm.gamma, smhm.lambda, smhm.scatter, smhm.mu, smhm.kappa);
  }
  //if (smf) phi_corr = 1.0 /(1.0+smhm.kappa);
  if (!smf) {
    
    // phi_corr *= pow(10,-sm_true)*smhm.ssfr_corr;
    // Apply an additional offset of 10^(kappa * exp(-(z-2)^2/2)) to SSFR.
    phi_corr *= pow(10,-sm_true+smhm.kappa * exp(-0.5 * (z - 2) *(z - 2)))*smhm.ssfr_corr;
    
    // // Apply a uniform offset of 10^(mu + kappa * exp(-(z-2)^2/2)) to SSFR.
    // phi_corr *= pow(10,-sm_true+smhm.mu+smhm.kappa * exp(-0.5 * (z - 2) *(z - 2)))*smhm.ssfr_corr;
    ih.corr = smhm.sfr_sm_corr;
  }
  else
  {
    phi_corr *= smhm.sm_completeness;
    ih.corr = 0;
  }

  if (!use_obs_psf) smhm.obs_scatter = 0;
  ih.gauss_inv_sigma = gen_inv_sigma(smhm.obs_scatter, smhm.scatter);
  // ih.gauss_inv_sigma = 10;
  // ih.kappa = smhm.kappa;
  // ih.passive_mass = smhm.passive_mass;
  ih.s_corr = scatter_corr;
  ih.sm = sm_true;
  ih.scatter = smhm.scatter;
  ih.obs_scatter = smhm.obs_scatter;
  ih.n = ih.n2 = -1;
  if (smf) 
  {
    ih.list = steps[ih.step].smf;
    ih.sfrac_sm = steps[ih.step].sfrac_sm;
    //if (z <= 0.2) {
      ih.n = ih.step;
      ih.ok = steps[ih.step].smf_ok;
      //}
    if (ih.f > -1) {
      ih.list2 = steps[ih.step+1].smf;
      // ih.sfrac_sm2 = steps[ih.step+1].sfrac_sm;
      //if (z <= 0.2) {
  ih.n2 = ih.step+1;
  ih.ok2 = steps[ih.step+1].smf_ok;
  // if (output_smf == 0)
  // {
  //   for (int j = SM_EXTRA; j < SM_BINS+SM_EXTRA; j++)
  //   {
  //     double sm_tmp = (SM_MIN + (j * 1.0 - SM_EXTRA) / SM_BPDEX);
  //     printf("sm=%f, smf1[%d]=%f, smf2[%d]=%f\n", sm_tmp, j, log10(ih.list[j]), j, log10(ih.list2[j]));
  //   }
  // }
  // output_smf = 1;
  //}
    }
  } else {
    ih.list = steps[ih.step].sfr_sm;
    ih.sfrac_sm = steps[ih.step].sfrac_sm;
    ih.sfr_sm_sf = steps[ih.step].sfr_sm_sf;
    ih.sfr_sm_q = steps[ih.step].sfr_sm_q;
    if (ih.f > -1) 
    {
      ih.list2 = steps[ih.step+1].sfr_sm;
      //ih.sfrac_sm2 = steps[ih.step].sfrac_sm2;
      ih.sfr_sm_sf2 = steps[ih.step+1].sfr_sm_sf;
      ih.sfr_sm_q2 = steps[ih.step+1].sfr_sm_q;
    }
  }
  sm_max = sm_true+GAUSS_CACHE_MAX/ih.gauss_inv_sigma;
  if (sm_max > steps[ih.step].smhm.sm_max)
    sm_max = steps[ih.step].smhm.sm_max;
  //if (sm_true >= sm_max) return 0;
  // fprintf(stderr, "sm_max=%f\n", sm_max);
  double result;
  // smhm.scatter = 0;
  // smhm.obs_scatter = 0;
  if (!smhm.scatter && (!use_obs_psf || !smhm.obs_scatter))
  {
    if (smf == 2) //quenched fractions
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
  else 
  {
    double precision = PHI_INTEGRAL_PRECISION;
     if (z <= 0.2) precision*=0.001;
     //fprintf(stderr, "sm=%f, a (sm_true) =%f, b (sm_max)=%f, phi_corr=%f\n", sm, sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma, sm_max, phi_corr);
    if (smf == 2)
    {
      // double epsilon_q = step_integral_helper_Q(0.5 * (sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma + sm_max), (void *)&ih);
      // double epsilon_sf = step_integral_helper_SF(0.5 * (sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma + sm_max), (void *)&ih);

      // epsilon_q *= (sm_max - (sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma)) * PHI_INTEGRAL_PRECISION;
      // epsilon_sf *= (sm_max - (sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma)) * PHI_INTEGRAL_PRECISION;

      double smf_q = (phi_corr * adaptiveGauss(step_integral_helper_Q, (void *)&ih,
          sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
            sm_max,
               precision,omp_get_thread_num()));
      double smf_sf = (phi_corr * adaptiveGauss(step_integral_helper_SF, (void *)&ih,
          sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
            sm_max,
               precision,omp_get_thread_num()));

      // double smf_q = (phi_corr * adaptiveSimpsons(step_integral_helper_Q, (void *)&ih,
      //     sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
      //       sm_max,
      //          epsilon_q,10));
      // double smf_sf = (phi_corr * adaptiveSimpsons(step_integral_helper_SF, (void *)&ih,
      //     sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma,
      //       sm_max,
      //          epsilon_sf,10));



      result = smf_q / (smf_q + smf_sf);
      // fprintf(stderr, "Interested SM: %.6f, Integral range: [%.6f, %.6f]; smf_q: %.6e, smf_sf: %.6e\n",
      //         sm, sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma, sm_max, smf_q, smf_sf);
    }
    else if (smf == 1)
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
      
     //result = 0; 
     // double val_max = -999;
     // int kk = 0;
     // for (double sm_tmp = sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma; sm_tmp < sm_max; sm_tmp+= 0.05)
     // {
     //  double val_tmp = log10(step_integral_helper(sm_tmp, (void *)&ih));
     //  if (val_tmp >= val_max) val_max = val_tmp;
     //  fprintf(stderr, "%f %f %f\n", sm_true, sm_tmp, val_tmp);
     //  // fprintf(stderr, "sm=%f, sm_tmp=%f, log(val)=%f\n", sm_true, sm_tmp, val_tmp);
     // }
     // printf("sm=%f, val_max=%f\n", sm, val_max);
  }
  
  return result;
}


double evaluate_from_step_uv(double z, double uv) {
  double phi_corr=1;
  struct step_integral_helper_data ih;
  calc_step_at_z(z, &(ih.step), &(ih.f));
  struct smf smhm = steps[ih.step].smhm;
  // double scatter_corr = smhm.scatter*smhm.scatter*(0.5*M_LN10);

  // if (!use_obs_psf) smhm.obs_scatter = 0;
  // ih.gauss_inv_sigma = gen_inv_sigma(smhm.obs_scatter, smhm.scatter);
  // ih.gauss_inv_sigma = 10;
  // ih.kappa = smhm.kappa;
  // ih.passive_mass = smhm.passive_mass;
  //ih.s_corr = scatter_corr;
  ih.sm = uv;
  //ih.scatter = smhm.scatter;
  //ih.obs_scatter = smhm.obs_scatter;
  ih.n = ih.n2 = -1;
  
  ih.list = steps[ih.step].uvlf;
  ih.sfrac_sm = steps[ih.step].std_uvlf;
  //if (z <= 0.2) {
  ih.n = ih.step;
  ih.ok = steps[ih.step].uvlf_ok;
  
  //for (int j=0; j<UV_BINS-1; j++) {
  //  double uv_tmp = UV_MIN + (double)j*UV_INV_BPMAG;   
  //  fprintf(stderr, "uv=%f, uvlf[%d]=%e, uvlf_next[%d]=%e, std_uvlf[%d]=%f\n", uv_tmp, j, ih.list[j+UV_EXTRA], j, ih.list2[j+UV_EXTRA], j, ih.sfrac_sm[j+UV_EXTRA]);
//}  

  //}
  if (ih.f > -1) 
  {
    ih.list2 = steps[ih.step+1].uvlf;
    ih.n2 = ih.step+1;
    ih.ok2 = steps[ih.step+1].uvlf_ok;
  }

  //for (int j=0; j<UV_BINS-1; j++) {
    //double uv_tmp = UV_MIN + (double)j*UV_INV_BPMAG;
    //fprintf(stderr, "uv=%f, uvlf[%d]=%e, uvlf_ok[%d]=%d, uvlf_next[%d]=%e, uvlf_next_ok[%d]=%d, std_uvlf[%d]=%f\n", uv_tmp, j, ih.list[j+UV_EXTRA], j, ih.ok[j+UV_EXTRA], j, ih.list2[j+UV_EXTRA], j, ih.ok2[j+UV_EXTRA], j, ih.sfrac_sm[j+UV_EXTRA]);
//}
  
  double uv_max = uv+UV_MAG_OFFSET;
  if (uv_max > steps[ih.step].smhm.uv_max)
    uv_max = steps[ih.step].smhm.uv_max;
  double uv_min = uv - UV_MAG_OFFSET;
  if (uv_min < steps[ih.step].smhm.uv_min)
    uv_min = steps[ih.step].smhm.uv_min;
  //if (sm_true >= sm_max) return 0;
  // fprintf(stderr, "sm_max=%f\n", sm_max);
  double result;
  double precision = PHI_INTEGRAL_PRECISION;
   if (z <= 0.2) precision*=0.001;
   //fprintf(stderr, "sm=%f, a (sm_true) =%f, b (sm_max)=%f, phi_corr=%f\n", sm, sm_true + GAUSS_CACHE_MIN/ih.gauss_inv_sigma, sm_max, phi_corr);
  

    result = (phi_corr * adaptiveGauss(step_integral_helper_uv, (void *)&ih,
        uv_min,
          uv_max,
             precision,omp_get_thread_num())); 
      

  return result;
}

double calc_ssfr(double sm, double z) {
  double ssfr = evaluate_from_step(z, sm, 0);
  //fprintf(stderr, "SSFR: %e (sm: %f; z:%f)\n", ssfr, sm, z);
  if (ssfr < 1e-15) ssfr = 1e-15;
  return ssfr;
}

double _quasar_lf_helper(int64_t mb, void *extra_info) {
  //Mi(z=2) = -5.26 - 2.5 log_10(n M_BH)
  //log10(n) = -(Mi + 5.26)/2.5 - log_10(M_BH)
  struct quasar_info *qi = extra_info;
  if (mb < 0) { mb = 0; }
  if (mb >= M_BINS) { mb = M_BINS-1; }
  double mbh = steps[qi->step].log_bh_mass[mb];
  double eta = -(qi->l + 5.26)/2.5 - mbh;
  double nd = steps[qi->step].t[mb];
  //double alpha = steps[qi->step].smhm.bh_alpha;
  double eta_0 = steps[qi->step].bh_eta[mb];
  double eta_frac = eta - eta_0;
  //printf("%f %f %f %f %e %f %f %f\n", qi->l, m, mbh, eta, nd, alpha, eta_0, eta_frac);
  if (eta_frac < steps[qi->step].ledd_min[mb] || eta_frac > steps[qi->step].ledd_max[mb]) return 0; 
  //double prob = exp10fc(eta_frac*alpha)*exp(-exp10fc(eta_frac));

  double bher_f = (eta_frac-steps[qi->step].ledd_min[mb])*steps[qi->step].ledd_bpdex[mb];
  int64_t bher_b = bher_f;
  bher_f -= bher_b;

  double p1 = steps[qi->step].bher_dist_full[mb*BHER_BINS+bher_b];
  double p2 = steps[qi->step].bher_dist_full[mb*BHER_BINS+bher_b+1];
  if (bher_b >= BHER_BINS-1) p2 = p1;
  // double f_CTK = steps[qi->step].f_CTK[mb];
  double prob = p1 + bher_f*(p2-p1);
  // mass-dependent modulation of the duty cycle.
  // double f_mass = pow(exp10((M_MIN + (mb + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[qi->step].smhm.dc_beta);
  // f_mass = f_mass < 1? f_mass : 1;
  double f_mass = exp((log10(steps[qi->step].bh_mass_avg[mb]) - steps[qi->step].smhm.dc_mbh) / steps[qi->step].smhm.dc_mbh_w);
  f_mass = f_mass / (1 + f_mass);
  double dc = steps[qi->step].smhm.bh_duty * f_mass;
  if (dc < 1e-4) dc = 1e-4;
  prob *= dc;
  // prob *= (1 - f_CTK);
  return prob*nd;
}

double _quasar_lf_helper_ctn(int64_t mb, void *extra_info) {
  //Mi(z=2) = -5.26 - 2.5 log_10(n M_BH)
  //log10(n) = -(Mi + 5.26)/2.5 - log_10(M_BH)
  struct quasar_info *qi = extra_info;
  if (mb < 0) { mb = 0; }
  if (mb >= M_BINS) { mb = M_BINS-1; }
  double mbh = steps[qi->step].log_bh_mass[mb];
  double eta = -(qi->l + 5.26)/2.5 - mbh;
  double nd = steps[qi->step].t[mb];
  //double alpha = steps[qi->step].smhm.bh_alpha;
  double eta_0 = steps[qi->step].bh_eta[mb];
  double eta_frac = eta - eta_0;
  //printf("%f %f %f %f %e %f %f %f\n", qi->l, m, mbh, eta, nd, alpha, eta_0, eta_frac);
  if (eta_frac < steps[qi->step].ledd_min[mb] || eta_frac > steps[qi->step].ledd_max[mb]) return 0; 
  //double prob = exp10fc(eta_frac*alpha)*exp(-exp10fc(eta_frac));

  double bher_f = (eta_frac-steps[qi->step].ledd_min[mb])*steps[qi->step].ledd_bpdex[mb];
  int64_t bher_b = bher_f;
  bher_f -= bher_b;

  double p1 = steps[qi->step].bher_dist_full[mb*BHER_BINS+bher_b];
  double p2 = steps[qi->step].bher_dist_full[mb*BHER_BINS+bher_b+1];
  if (bher_b >= BHER_BINS-1) p2 = p1;
  // double f_CTK = steps[qi->step].f_CTK[mb];
  double prob = p1 + bher_f*(p2-p1);
  // mass-dependent modulation of the duty cycle.
  // double f_mass = pow(exp10((M_MIN + (mb + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[qi->step].smhm.dc_beta);
  // f_mass = f_mass < 1? f_mass : 1;
  double f_mass = exp((log10(steps[qi->step].bh_mass_avg[mb]) - steps[qi->step].smhm.dc_mbh) / steps[qi->step].smhm.dc_mbh_w);
  f_mass = f_mass / (1 + f_mass);
  double dc = steps[qi->step].smhm.bh_duty * f_mass;
  if (dc < 1e-4) dc = 1e-4;
  prob *= dc;
  // prob *= (1 - f_CTK);
  return prob*nd;
}

double quasar_lf_helper(double m, void *extra_info) {
  //Mi(z=2) = -5.26 - 2.5 log_10(n M_BH)
  //log10(n) = -(Mi + 5.26)/2.5 - log_10(M_BH)
  struct quasar_info *qi = extra_info;
  double mf = BPDEX*(m - M_MIN)-0.5;
  int64_t mb = mf;
  mf -= mb;
  if (mb < 0) { mb = 0; mf = 0; }
  if (mb >= M_BINS) { mb = M_BINS-2; mf = 1; }
  double mbh1 = steps[qi->step].log_bh_mass[mb];
  double mbh2 = steps[qi->step].log_bh_mass[mb+1];
  double mbh = mbh1 + mf*(mbh2-mbh1);
  double eta = -(qi->l + 5.26)/2.5 - mbh;
  double nd = exp10fc(mf_cache(steps[qi->step].scale, m));
  //double alpha = steps[qi->step].smhm.bh_alpha;
  double eta1 = steps[qi->step].bh_eta[mb];
  double eta2 = steps[qi->step].bh_eta[mb+1];
  double eta_0 = eta1 + mf*(eta2-eta1);
  double eta_frac = eta - eta_0;
  //printf("%f %f %f %f %e %f %f %f\n", qi->l, m, mbh, eta, nd, alpha, eta_0, eta_frac);
  if (eta_frac < steps[qi->step].ledd_min[mb] || eta_frac < steps[qi->step].ledd_min[mb+1] ||
      eta_frac > steps[qi->step].ledd_max[mb] || eta_frac > steps[qi->step].ledd_max[mb+1])
    return 0;
  //double prob = exp10fc(eta_frac*alpha)*exp(-exp10fc(eta_frac));

  double bher_f1 = (eta_frac-steps[qi->step].ledd_min[mb])*steps[qi->step].ledd_bpdex[mb];
  int64_t bher_b1 = bher_f1;
  bher_f1 -= bher_b1;

  double bher_f2 = (eta_frac-steps[qi->step].ledd_min[mb+1])*steps[qi->step].ledd_bpdex[mb+1];
  int64_t bher_b2 = bher_f2;
  bher_f2 -= bher_b2;

  /*
  double p1 = steps[qi->step].bher_dist[bher_b];
  double p2 = steps[qi->step].bher_dist[bher_b+1];
  if (bher_b >= BHER_BINS-1) p2 = p1;
  double prob = p1 + bher_f*(p2-p1);
  */

  double p1m1 = steps[qi->step].bher_dist_full[mb*BHER_BINS+bher_b1];
  double p2m1 = steps[qi->step].bher_dist_full[mb*BHER_BINS+bher_b1+1];
  if (bher_b1 >= BHER_BINS-1) p2m1 = p1m1;
  double pm1 = p1m1 + bher_f1*(p2m1-p1m1);
  double p1m2 = steps[qi->step].bher_dist_full[(mb+1)*BHER_BINS+bher_b2];
  double p2m2 = steps[qi->step].bher_dist_full[(mb+1)*BHER_BINS+bher_b2+1];
  if (bher_b2 >= BHER_BINS-1) p2m2 = p1m2;
  double pm2 = p1m2 + bher_f2*(p2m2-p1m2);
  double prob = pm1 + mf*(pm2-pm1);

  return prob*nd;
}

// double calc_bhmf(double m, double z) {
//   int64_t step;
//   double f;
//   calc_step_at_z(z, &step, &f);
//   int64_t i;
//   if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
//   double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.gamma);
//   double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
//   double s3 = s1*s1+s2*s2;
//   if (!s3) return -1;
//   double nd = 0;

//   double bhm[M_BINS]={0}, t[M_BINS]={0};
//   for (i=0; i<M_BINS; i++) {
//     bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
//     double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
//     t[i] = ndm;
//   }

//   double tw = 0;
//   int64_t s;
//   for (s=0; s<=(M_BINS-1)*5; s++) {
//     i = s/5;
//     f = (s%5)/5.0;
//     if (i==(M_BINS-1)) { i--; f=1; }
//     double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
//     double tnd = t[i] + f*(t[i+1]-t[i]);
//     double dm = tbhm - m;
//     double weight = exp(-0.5*dm*dm/s3);
//     tw += weight;
//     nd += weight*tnd;
//   }
//   nd /= tw*INV_BPDEX;
//   if (nd > 1e-15) return nd;
//   return 1e-15;
// }



double calc_bhmf(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  //double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.gamma);
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0};
  for (i=0; i<M_BINS; i++) {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
  }
  // nd /= tw*INV_BPDEX;
  nd /= 5.0;
  if (nd > 1e-15) return nd;
  return 1e-15;
}


double calc_bhmf_typeI(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }

  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  double bhm_f = (m - MBH_MIN) * mbh_bpdex;
  int64_t bhm_b = bhm_f;
  bhm_f -= bhm_b;

  if (bhm_b >= MBH_BINS - 1) {bhm_b = MBH_BINS - 2; bhm_f = 1;}
  double bhmf_tot1 = steps[step].bhmf[bhm_b] + bhm_f * (steps[step].bhmf[bhm_b+1] - steps[step].bhmf[bhm_b]);
  double bhmf_tot2 = steps[step+1].bhmf[bhm_b] + bhm_f * (steps[step+1].bhmf[bhm_b+1] - steps[step+1].bhmf[bhm_b]);
  double bhmf_tot = bhmf_tot1 + f * (bhmf_tot2 - bhmf_tot1);


  double f_mass = exp((m - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
  f_mass = f_mass / (1 + f_mass);
  double dc = steps[step].smhm.bh_duty * f_mass;
  if (dc < 1e-4) dc = 1e-4;

  double corr_obs = 0;
  double tw = 0;
  for (i=0; i<LBOL_BINS; i++)
  {
    double lbol = LBOL_MIN + (i + 0.5) * LBOL_INV_BPDEX;
    double lx = find_Lx_at_Lbol(lbol);
    double f_obs = F_obs(lx);
    corr_obs += (1 - f_obs) * steps[step].lum_dist_full[i];
    tw += steps[step].lum_dist_full[i];
  }  
  if (tw > 0) corr_obs /= tw;
  return bhmf_tot * corr_obs * dc;

}


double calc_active_bhmf(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  // pfile = fopen('./mcmc_runs/ABHMF_check.txt', 'w');
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*steps[step].smhm.scatter*steps[step].smhm.bh_gamma + f*(steps[step+1].smhm.scatter)*steps[step+1].smhm.bh_gamma;
  // double s1 = 0.3;
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double abhmf_shift = (1.0-f)*steps[step].smhm.abhmf_shift + f*(steps[step+1].smhm.abhmf_shift);
  if (z > 0.5) abhmf_shift = 0;
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhm[M_BINS]={0}, t[M_BINS]={0};
  for (i=0; i<M_BINS; i++) {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    double nactive = (1.0-f)*steps[step].f_active[i] + f*steps[step+1].f_active[i];
    //double f_CTK = steps[step].f_CTK[i];
    //t[i] = ndm*nactive * (1 - f_CTK);
    t[i] = ndm*nactive;
    //fprintf(stderr, "%e %e %e %e %e %e\n", m, M_MIN + (i + 0.5) * INV_BPDEX, bhm[i], s3, ndm, nactive);
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double dm = tbhm - (m + abhmf_shift);
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
  }
  // nd /= tw*INV_BPDEX;
  nd /= 5.0;
  if (nd > 1e-15) return nd;
  return 1e-15;
}

// calculate the average BHAR as a function of Mbh and z.
double calc_bhar_mbh(double m, double z) {
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
  double bhar_nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0};
  double bhm_avg[M_BINS] = {0};
  double bhar[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    if (!isfinite(tbhm_avg)) continue;
    //fprintf(stderr, "tbhm=%f, tbhm_avg=%f\n", tbhm, tbhm_avg);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    tbhar *= exp10fc(m - tbhm_avg); //constant BHER distribution!
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  bhar_nd /= nd;
  if (bhar_nd > 1e-15) return bhar_nd;
  return 1e-15;
}

// calculate the average BHAR as a function of Mstar and z.
double calc_bhar_mstar(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step].smhm.scatter);
  double mu = (1.0-f)*(steps[step].smhm.mu) + f*(steps[step].smhm.mu);
  // double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // double s3 = s1*s1+s2*s2;
  double s3 = s1*s1;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;


  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);

  double sm[M_BINS]={0}, t[M_BINS]={0}, mbh_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    mbh_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tmbh_avg = mbh_avg[i] + f*(mbh_avg[i+1]-mbh_avg[i]);
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    tbhar *= exp10fc(m - tmbh_avg);
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  bhar_nd /= nd;
  if (bhar_nd > 1e-15) return bhar_nd;
  return 1e-15;
}


// calculate the average BHAR/SFR ratio as a function of Mbh and z.
// This function is actually more complicated than the calculation of
// BHAR, BHMR, and BHERs, because there's a correlation between the 
// SFR and the average stellar mass. Consequently we need to trace
// the galaxy mass distribution in each halo mass, ***given the
// BH mass*** that we're looking for.
double calc_bhar_sfr_mbh(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
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
  for (i=0; i<M_BINS; i++) {
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
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg)) continue;
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double sfr_fine_bin = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(m - tbhm_avg);

    // double tsm_avg = 0;
    double tsfr = 0;
    double tnorm = 0;

    for (double sm_tmp=3; sm_tmp<=13; sm_tmp+=0.05)
    {
      double bm_tmp = bulge_mass(sm_tmp + mu, 1/(1+z));
      double dbhm_tmp = (m - calc_bh_at_bm(bm_tmp, steps[step].smhm)) / s2;
      double dsm_tmp = (sm_tmp - tsm) / scatter;
      double sfr_tmp = sfr_fine_bin * pow(exp10(sm_tmp) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);
      double prob = scatter_tot / (s2 * scatter) * exp(-0.5*dbhm_tmp*dbhm_tmp - 0.5*dsm_tmp*dsm_tmp);
      tnorm += prob;
      tsfr += sfr_tmp * prob;
      //fprintf(stderr, "i=%d, f=%f, sm_tmp=%f, bm_tmp=%f, dbhm_tmp=%f, dsm_tmp=%f, sfr_tmp=%e, prob=%e, tnorm=%e, tsfr=%e\n", i, f, sm_tmp, bm_tmp, dbhm_tmp, dsm_tmp, sfr_tmp, prob, tnorm, tsfr);
      // tsm_avg += exp10(sm_tmp) * prob;
    }
    // tsm_avg /= tnorm;
    tsfr /= tnorm;


    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
    sfr_nd += weight*tnd*tsfr;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  // bhar_nd /= nd;
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
double calc_bhar_sfr_mstar(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
  // double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // double scatter = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
  // double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double s3 = s1 * s1;
  // double s3 = s1*s1+s2*s2;
  // double scatter_tot = sqrt(s3);
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhar_nd = 0;
  double sfr_nd = 0;

  double mbulge = bulge_mass(m + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(mbulge, steps[step].smhm);

  double sm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (s1 * M_LN10 * s1 * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg)) continue;
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(mbh_med - tbhm_avg);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);

    // // double tsm_avg = 0;
    // double tsfr = 0;
    // double tnorm = 0;
    // for (double sm_tmp=3, sm_tmp<=13, sm_tmp+=0.05)
    // {
    //   double bm_tmp = bulge_mass(sm_tmp + mu, 1/(1+z));
    //   double dbhm_tmp = (m - calc_bh_at_bm(bm_tmp, steps[step].smhm)) / s2;
    //   double dsm_tmp = (sm_tmp - tsm) / scatter;
    //   double sfr_tmp = tsfr * pow(sm_tmp / (tsm * exp(0.5 * (scatter * M_LN10 * scatter * M_LN10))), sfr_sm_corr);
    //   double prob = scatter_tot / (s2 * scatter) * exp(-0.5*dbhm_tmp*dbhm_tmp - 0.5*dsm_tmp*dsm_tmp);
    //   tnorm += prob;
    //   tsfr += sfr_tmp * prob;
    //   // tsm_avg += exp10(sm_tmp) * prob;
    // }
    // // tsm_avg /= tnorm;
    // tsfr /= norm;


    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
    sfr_nd += weight*tnd*tsfr;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  // bhar_nd /= nd;
  double bhar_sfr = bhar_nd / sfr_nd;
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}


double calc_sfr_mstar(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
  // double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // double scatter = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
  // double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1 * s1;
  // double s3 = s1*s1+s2*s2;
  // double scatter_tot = sqrt(s3);
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  // double bhar_nd = 0;
  double sfr_nd = 0;

  double sm[M_BINS]={0}, t[M_BINS]={0};
  // double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    // bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    // bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (s1 * M_LN10 * s1 * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    // double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    // double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);



    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    // bhar_nd += weight*tnd*tbhar;
    sfr_nd += weight*tnd*tsfr;
  }
  sfr_nd /= nd;
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  // bhar_nd /= nd;
  
  if (sfr_nd > 1e-15) return sfr_nd;
  return 1e-15;
}




// calculate the average SBHAR/SSFR ratio as a function of Mbh and z.
// This function is actually more complicated than the calculation of
// BHAR, BHMR, and BHERs, because there's a correlation between the 
// SFR and the average stellar mass. Consequently we need to trace
// the galaxy mass distribution in each halo mass, ***given the
// BH mass*** that we're looking for.
double calc_sbhar_ssfr_mbh(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
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
  for (i=0; i<M_BINS; i++) {
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
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    if (sm[i+1] == 0) f=0;
    if (sm[i] == 0) continue;
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
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
      //fprintf(stderr, "i=%d, f=%f, tsm=%f, sm_tmp=%f, bm_tmp=%f, dbhm_tmp=%f, dsm_tmp=%f, sfr_tmp=%e, prob=%e, tnorm=%e, tsfr=%e\n", i, f, tsm, sm_tmp, bm_tmp, dbhm_tmp, dsm_tmp, sfr_tmp, prob, tnorm, tsfr);
    }
    tsm_avg /= tnorm;
    tsfr /= tnorm;
    //fprintf(stderr, "i=%d, f=%f, tsm for the bin=%f, np.log10(avg sm at fixed mbh / sm_corr)=%f\n", i, f, tsm, log10(tsm_avg / med_to_avg_star));

    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
    sfr_nd += weight*tnd*tsfr;
    sm_avg_nd += weight * tnd * tsm_avg;
  }
  //nd /= 5.0;
  bhar_nd /= nd; //divide by the number density
  sfr_nd /= nd;
  sm_avg_nd /= nd;
  double bhar_sfr = bhar_nd / sfr_nd / (exp10(m)) * sm_avg_nd;
  double bm_scaling = (m - steps[step].smhm.bh_beta) / steps[step].smhm.bh_gamma + 11;
  // double sm_scaling = stellar_mass(bm_scaling, 1 / (1 + z)) - mu;
  // fprintf(stderr, "mbh=%e, z=%f,  mstar_from_scaling=%f, mstar_from_integral=%f, dMstar=%f, log10(bhar/sfr)=%f\n", m, z, sm_scaling, log10(sm_avg_nd/med_to_avg_star), sm_scaling - log10(sm_avg_nd/med_to_avg_star), log10(bhar_nd/sfr_nd)); 
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}




// calculate the average SBHAR/SSFR ratio as a function of Mstar and z.
// This function is actually more complicated than the calculation of
// BHAR, BHMR, and BHERs, because there's a correlation between the 
// SFR and the average stellar mass. Consequently we need to trace
// the galaxy mass distribution in each halo mass, ***given the
// BH mass*** that we're looking for.
double calc_sbhar_ssfr_mstar(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.bh_gamma) + f*(steps[step+1].smhm.scatter * steps[step+1].smhm.bh_gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  double scatter = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step+1].smhm.scatter);
  double mu = (1.0-f)*steps[step].smhm.mu + f*(steps[step+1].smhm.mu);
  double sfr_sm_corr = (1.0-f)*steps[step].smhm.sfr_sm_corr + f*(steps[step+1].smhm.sfr_sm_corr);
  double s3 = s1 * s1;
  // double s3 = s1*s1+s2*s2;
  // double scatter_tot = sqrt(s3);
  if (!s3) 
  {
    fprintf(stderr, "s3=%f, s1=%f, scatter=%f, gamma=\n");
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

  double sm[M_BINS]={0}, t[M_BINS]={0}, bhm_avg[M_BINS]={0};
  double bhar[M_BINS] = {0};
  double sfr[M_BINS] = {0};
  double bhm[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    bhm_avg[i] = (1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]);
    bhar[i] = (1.0-f)*steps[step].bh_acc_rate[i] + f*steps[step+1].bh_acc_rate[i];
    sfr[i] = (1.0-f)*steps[step].sfr[i] + f*steps[step+1].sfr[i];
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    // sm_avg[i] = (1.0-f)*steps[step].sm_avg[i] + f*steps[step+1].sm_avg[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  double med_to_avg_star = exp(0.5 * (scatter * M_LN10 * scatter * M_LN10));
  double med_to_avg_bh = exp(0.5 * ((s2 * s2) * M_LN10 * M_LN10));
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm_avg = bhm_avg[i] + f*(bhm_avg[i+1]-bhm_avg[i]);
    if (!isfinite(tbhm_avg)) continue;
    double tbm = bulge_mass(m + mu, 1/(1+z));
    double tbhm = calc_bh_at_bm(tbm, steps[step].smhm);
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    // double tsm_avg = sm_avg[i] + f*(sm_avg[i+1]-sm_avg[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhar = bhar[i] + f*(bhar[i+1]-bhar[i]);
    double tsfr = sfr[i] + f*(sfr[i+1]-sfr[i]);

    tbhar *= exp10fc(mbh_med - tbhm_avg);

    tsfr *= pow(exp10(m) / (exp10(tsm) * med_to_avg_star), sfr_sm_corr);

    // // double tsm_avg = 0;
    // double tsfr = 0;
    // double tnorm = 0;
    // for (double sm_tmp=3, sm_tmp<=13, sm_tmp+=0.05)
    // {
    //   double bm_tmp = bulge_mass(sm_tmp + mu, 1/(1+z));
    //   double dbhm_tmp = (m - calc_bh_at_bm(bm_tmp, steps[step].smhm)) / s2;
    //   double dsm_tmp = (sm_tmp - tsm) / scatter;
    //   double sfr_tmp = tsfr * pow(sm_tmp / (tsm * exp(0.5 * (scatter * M_LN10 * scatter * M_LN10))), sfr_sm_corr);
    //   double prob = scatter_tot / (s2 * scatter) * exp(-0.5*dbhm_tmp*dbhm_tmp - 0.5*dsm_tmp*dsm_tmp);
    //   tnorm += prob;
    //   tsfr += sfr_tmp * prob;
    //   // tsm_avg += exp10(sm_tmp) * prob;
    // }
    // // tsm_avg /= tnorm;
    // tsfr /= norm;


    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhar_nd += weight*tnd*tbhar;
    sfr_nd += weight*tnd*tsfr;
    // sm_avg_nd += weight*tnd*exp10(tsm) * med_to_avg_star;
    bhm_avg_nd += weight*tnd*exp10(tbhm) * med_to_avg_bh;
  }
  // nd /= tw*INV_BPDEX;
  //nd /= 5.0;
  bhar_nd /= nd;
  sfr_nd /= nd;
  // sm_avg_nd /= nd;
  bhm_avg_nd /= nd;
  //fprintf(stderr, "m=%f, z=%f, bhar=%e, sfr=%e, bhm=%e, sm_corr=%e, bh_corr=%e\n", m, z, bhar_nd, sfr_nd, bhm_avg_nd, med_to_avg_star, med_to_avg_bh);
  double bhar_sfr = bhar_nd / sfr_nd / (bhm_avg_nd / (exp10(m)));
  double tbm = bulge_mass(m + mu, 1/(1+z));
  double tbhm = calc_bh_at_bm(tbm, steps[step].smhm);
  fprintf(stderr, "mstar=%f, z=%f, mbh_from_scaling=%f, mbh_from_integral=%f, dmbh=%f, log10(bhar/sfr)=%f\n", m, z, tbhm, log10(bhm_avg_nd/med_to_avg_bh), tbhm - log10(bhm_avg_nd/med_to_avg_bh), log10(bhar_nd/sfr_nd));
  if (bhar_sfr > 1e-15) return bhar_sfr;
  return 1e-15;
}



// calculate the average BHER as a function of Mbh and z.
double calc_bher_mbh(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter * steps[step].smhm.gamma) + f*(steps[step].smhm.scatter * steps[step].smhm.gamma);
  double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // double eff_rad = (1.0-f)*steps[step].smhm.bh_efficiency_rad + f*(steps[step+1].smhm.bh_efficiency_rad);
  double s3 = s1*s1+s2*s2;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bher_nd = 0;

  double bhm[M_BINS]={0}, t[M_BINS]={0};
  double bher[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];

    // steps[n].bh_eta[i] = 
    // log10(//schechter_inv_avg(steps[n].smhm.bh_alpha, BHER_MIN, BHER_EFF_MAX)*
    // bhar_tmp/steps[n].bh_mass_avg[i]*4.5e8*(steps[n].smhm.bh_efficiency_rad));


    bher[i] = 4.5e8 * ((1.0-f)*(steps[step].bh_acc_rate[i] / steps[step].bh_mass_avg[i] * steps[step].smhm.bh_efficiency_rad) + 
                  f*steps[step+1].bh_acc_rate[i] / steps[step+1].bh_mass_avg[i] * steps[step+1].smhm.bh_efficiency_rad);
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbher = bher[i] + f*(bher[i+1]-bher[i]);
    if (!isfinite(tbher)) continue;
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bher_nd += weight*tnd*tbher;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  //fprintf(stderr, "m=%f, bher_nd=%e, nd=%e, bher=%e\n", m, bher_nd, nd, bher_nd / nd);
  bher_nd /= nd;
  
  if (bher_nd > 1e-15) return bher_nd;
  return 1e-15;
}

// calculate the average BHER as a function of Mstar and z.
double calc_bher_mstar(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step].smhm.scatter);
  // double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // double eff_rad = (1.0-f)*steps[step].smhm.bh_efficiency_rad + f*(steps[step+1].smhm.bh_efficiency_rad);
  // double s3 = s1*s1+s2*s2;
  double s3 = s1*s1;
  //fprintf(stderr, "scatter=%f\n", s3);
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bher_nd = 0;

  double sm[M_BINS]={0}, t[M_BINS]={0};
  double bher[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];

    // steps[n].bh_eta[i] = 
    // log10(//schechter_inv_avg(steps[n].smhm.bh_alpha, BHER_MIN, BHER_EFF_MAX)*
    // bhar_tmp/steps[n].bh_mass_avg[i]*4.5e8*(steps[n].smhm.bh_efficiency_rad));


    bher[i] = 4.5e8 * ((1.0-f)*(steps[step].bh_acc_rate[i] / steps[step].bh_mass_avg[i] * steps[step].smhm.bh_efficiency_rad) + 
                  f*steps[step+1].bh_acc_rate[i] / steps[step+1].bh_mass_avg[i] * steps[step+1].smhm.bh_efficiency_rad);
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbher = bher[i] + f*(bher[i+1]-bher[i]);
    if (!isfinite(tbher)) continue;
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bher_nd += weight*tnd*tbher;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  //fprintf(stderr, "m=%f, bher_nd=%e, nd=%e, bher=%e\n", m, bher_nd, nd, bher_nd / nd);
  bher_nd /= nd;
  if (bher_nd > 1e-15) return bher_nd;
  return 1e-15;
}



// calculate the average BHMR as a function of Mbh and z.
double calc_bhmr_mbh(double m, double z) {
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

  double bhm[M_BINS]={0}, t[M_BINS]={0};
  double bhmr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    bhm[i] = (1.0-f)*steps[step].log_bh_mass[i] + f*steps[step+1].log_bh_mass[i];
    bhmr[i] = (1.0-f)*steps[step].bh_merge_rate[i] + f*steps[step+1].bh_merge_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tbhm = bhm[i] + f*(bhm[i+1]-bhm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhmr = bhmr[i] + f*(bhmr[i+1]-bhmr[i]);
    double dm = tbhm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhmr_nd += weight*tnd*tbhmr;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  bhmr_nd /= nd;
  if (bhmr_nd > 1e-15) return bhmr_nd;
  return 1e-15;
}

// calculate the average BHMR as a function of Mstar and z.
double calc_bhmr_mstar(double m, double z) {
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  int64_t i;
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = (1.0-f)*(steps[step].smhm.scatter) + f*(steps[step].smhm.scatter);
  // double s2 = (1.0-f)*steps[step].smhm.bh_scatter + f*(steps[step+1].smhm.bh_scatter);
  // double s3 = s1*s1+s2*s2;
  double s3 = s1*s1;
  if (!s3) return 1e-15;
  double norm_gauss = 1 / sqrt(2 * M_PI * s3);
  double nd = 0;
  double bhmr_nd = 0;

  double sm[M_BINS]={0}, t[M_BINS]={0};
  double bhmr[M_BINS] = {0};
  for (i=0; i<M_BINS; i++) {
    sm[i] = (1.0-f)*steps[step].log_sm[i] + f*steps[step+1].log_sm[i];
    bhmr[i] = (1.0-f)*steps[step].bh_merge_rate[i] + f*steps[step+1].bh_merge_rate[i];
    double ndm = (1.0-f)*steps[step].t[i] + f*steps[step+1].t[i];
    t[i] = ndm;
  }

  // double tw = 0;
  int64_t s;
  for (s=0; s<=(M_BINS-1)*5; s++) {
    i = s/5;
    f = (s%5)/5.0;
    if (i==(M_BINS-1)) { i--; f=1; }
    double tsm = sm[i] + f*(sm[i+1]-sm[i]);
    double tnd = t[i] + f*(t[i+1]-t[i]);
    double tbhmr = bhmr[i] + f*(bhmr[i+1]-bhmr[i]);
    double dm = tsm - m;
    double weight = exp(-0.5*dm*dm/s3) * norm_gauss;
    // tw += weight;
    nd += weight*tnd;
    bhmr_nd += weight*tnd*tbhmr;
  }
  // nd /= tw*INV_BPDEX;
  // nd /= 5.0;
  bhmr_nd /= nd;
  if (bhmr_nd > 1e-15) return bhmr_nd;
  return 1e-15;
}


double cosmic_bh_density(double z, double thresh_ledd, double thresh_lbol, struct smf_fit *fit) {
  int64_t i, step;
  double f;
  calc_step_at_z(z, &step, &f);
  double mnd = 0;
  double alpha = (1.0-f)*steps[step].smhm.bh_alpha + f*(steps[step+1].smhm.bh_alpha);
  double delta = (1.0-f)*steps[step].smhm.bh_delta + f*(steps[step+1].smhm.bh_delta);
  double bh_eta_crit = (1.0-f)*steps[step].smhm.bh_eta_crit + f*(steps[step+1].smhm.bh_eta_crit);
  for (i=0; i<M_BINS; i++) {
    double dc = (1.0-f)*steps[step].smhm.bh_duty + f*(steps[step+1].smhm.bh_duty);
    // dc *= pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
    // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
    // f_mass = f_mass < 1? f_mass : 1;
    double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;
    double sn = dc/doublePL_frac_above_thresh(BHER_EFF_MIN, alpha, delta, 1.0, fit);
    double bhm = (1.0-f)*steps[step].bh_mass_avg[i] + f*(steps[step+1].bh_mass_avg[i]);
    //double bhm = exp10((1.0-f)*steps[step].log_bh_mass[i] + f*(steps[step+1].log_bh_mass[i]));
    //double bhm = exp10((1.0-f)*log10(steps[step].bh_mass_avg[i]) + f*log10(steps[step+1].bh_mass_avg[i]));
    double nd = (1.0-f)*steps[step].t[i] + f*(steps[step+1].t[i]);
    if (thresh_lbol) {
      double lbhm = (1.0-f)*steps[step].log_bh_mass[i] + f*(steps[step+1].log_bh_mass[i]);
      double lum_i = 72.5 - 2.5*(thresh_lbol - 7);
      thresh_ledd = -(lum_i + 5.26)/2.5 - lbhm;
      if (nonlinear_luminosity) {
  const double log_of_2 = M_LN2 / M_LN10;
  if (thresh_ledd > log_of_2) {
    double bher_norm = exp10(thresh_ledd - log_of_2);
    if (!isfinite(bher_norm)) continue;
    thresh_ledd = (bher_norm - 1.0 + M_LN2)/M_LN10;
  }
  else if (thresh_ledd < bh_eta_crit) {
    thresh_ledd = (thresh_ledd + bh_eta_crit)/2.0;
  }
      }
    }

    if (thresh_ledd) {
      double beta = (1.0-f)*steps[step].bh_eta[i] + f*(steps[step+1].bh_eta[i]);
      double t_ef = thresh_ledd - beta;
      double f_active = doublePL_frac_above_thresh(t_ef, alpha, delta, sn, fit);
      mnd += nd*bhm*f_active;
    } else {
      mnd += nd*bhm;
    }
  }
  return mnd;
}

double cosmic_bh_density_split(double z, double Mh_low, double Mh_high, struct smf_fit *fit) {
  int64_t i, step;
  double f;
  double mh;
  calc_step_at_z(z, &step, &f);
  double mnd = 0;
  double nd_tot = 0;
  // double alpha = (1.0-f)*steps[step].smhm.bh_alpha + f*(steps[step+1].smhm.bh_alpha);
  // double delta = (1.0-f)*steps[step].smhm.bh_delta + f*(steps[step+1].smhm.bh_delta);
  // double bh_eta_crit = (1.0-f)*steps[step].smhm.bh_eta_crit + f*(steps[step+1].smhm.bh_eta_crit);
  for (i=0; i<M_BINS; i++) {
    mh = M_MIN + (i + 0.5) * INV_BPDEX;
    if (mh < Mh_low || mh > Mh_high) continue;
    // double dc = (1.0-f)*steps[step].smhm.bh_duty + f*(steps[step+1].smhm.bh_duty);
    // dc *= pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
    // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
    // f_mass = f_mass < 1? f_mass : 1;
    // double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    // f_mass = f_mass / (1 + f_mass);
    // dc *= f_mass;
    // if (dc < 1e-4) dc = 1e-4;
    // double sn = dc/doublePL_frac_above_thresh(BHER_EFF_MIN, alpha, delta, 1.0, fit);
    //double bhm = (1.0-f)*steps[step].bh_mass_avg[i] + f*(steps[step+1].bh_mass_avg[i]);
    double bhm = exp10((1.0-f)*steps[step].log_bm[i] + f*(steps[step+1].log_bm[i]));
    double nd = (1.0-f)*steps[step].t[i] + f*(steps[step+1].t[i]);
    mnd += nd*bhm;
    nd_tot += nd;
  }
  mnd /= nd_tot;

  return mnd;
}

double calc_quasar_lf(double l, double z) {
  struct quasar_info qi;
  qi.l = l;
  double sf = 0;
  calc_step_at_z(z, &(qi.step), &sf);
  int64_t i=0;


  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  double ld = 0;
  for (i=0; i<M_BINS; i++)
    ld += _quasar_lf_helper(i, &qi);
  return norm*ld;
}

double calc_quasar_lf_new(double Mi, double z) {
  double lbol = 36 - 0.4 * Mi;
  int64_t step; double sf = 0;
  calc_step_at_z(z, &(step), &sf);
  int64_t i=0;
  double mbh_min = steps[step].bh_mass_min; double mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  if (lbol < LBOL_MIN || lbol > LBOL_MAX) return 0;
  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  double ld = 0;
  for (i=0; i<MBH_BINS; i++)
  {
    double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
    //fprintf(stderr, "lbol=%f, LBOL_MIN=%f, LBOL_BPDEX=%f\n", lbol, LBOL_MIN, LBOL_BPDEX);
    int64_t lbol_b = lbol_f;
    lbol_f -= lbol_b;
    double p1 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b];
    double p2 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b+1];
    ld += p1 + lbol_f * (p2 - p1);
    //fprintf(stderr, "step=%d, i=%d, Mbh=%f, p1=%e. p2=%e, lbol_f=%f, lbol_b=%d\n", step, i, mbh_min + (i + 0.5) * mbh_inv_bpdex, p1, p2, lbol_f, lbol_b);
  }
  return norm*ld;
}

double calc_quasar_lf_mbh(double Mi, double z, double mbh_low, double mbh_high) {
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

// To calculate the QLF in slices of Eddington ratios, we simply convert the eta limits
// to Mbh limits.
double calc_quasar_lf_eta(double Mi, double z, double eta_low, double eta_high) {
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

double calc_quasar_lf_ctn(double l, double z) {
  struct quasar_info qi;
  qi.l = l;

  // convert l back to Log10(erg/s) in order to calculate the Compton-thick fraction.
  l = (90 - l) / 2.5;
  double Lx = find_Lx_at_Lbol(l);
  double psi_ctn = psi(Lx, z);

  double sf = 0;
  calc_step_at_z(z, &(qi.step), &sf);
  int64_t i=0;

  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex

  double ld = 0;
  for (i=0; i<M_BINS; i++)
    ld += _quasar_lf_helper_ctn(i, &qi);
  return norm*ld / (1 + psi_ctn); //The Compton-thick correction comes from that we assume the fraction
  // of CTK objects is the same as CTN ones in each redshift and luminosity.
}

double calc_quasar_lf_split(double l, double z, double Mh_min, double Mh_max) {
  struct quasar_info qi;
  qi.l = l;
  double sf = 0;
  calc_step_at_z(z, &(qi.step), &sf);
  int64_t i=0;


  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex
  // norm *= steps[qi.step].smhm.bh_duty;

  /*  return norm * adaptiveGauss(quasar_lf_helper, (void *)&qi,
            m_min, m_max, PHI_INTEGRAL_PRECISION*0.001, 
            omp_get_thread_num());
  */

  double ld = 0;
  for (i=0; i<M_BINS; i++)
  {
    if (M_MIN + (i + 0.5) * INV_BPDEX < Mh_min || M_MIN + (i + 0.5) * INV_BPDEX > Mh_max) continue;
    ld += _quasar_lf_helper(i, &qi);
  }
    
  return norm*ld;
}

double calc_quasar_lf_split_sm(double l, double z, double sm_min, double sm_max) {
  struct quasar_info qi;
  qi.l = l;
  double sf = 0;
  calc_step_at_z(z, &(qi.step), &sf);
  int64_t i=0;


  double norm = 1.0;  //1.0/schechter_norm(steps[qi.step].smhm.bh_alpha, BHER_MIN, BHER_MAX);
  norm /= 2.5; //Per mag instead of per dex
  // norm *= steps[qi.step].smhm.bh_duty;

  /*  return norm * adaptiveGauss(quasar_lf_helper, (void *)&qi,
            m_min, m_max, PHI_INTEGRAL_PRECISION*0.001, 
            omp_get_thread_num());
  */

  double ld = 0;
  for (i=0; i<M_BINS; i++)
  {
    if (log10(steps[qi.step].sm_avg[i]) < sm_min || log10(steps[qi.step].sm_avg[i]) > sm_max) continue;
    ld += _quasar_lf_helper(i, &qi);
  }
    
  return norm*ld;
}

double _prob_of_ledd_nonlinear(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit) {
  double bher = ledd;
  double bher_norm = 1;
  const double log_of_2 = M_LN2 / M_LN10;
  if (bher > log_of_2) {
    bher_norm = exp10(ledd - log_of_2);
    if (!isfinite(bher_norm)) return 0;
    bher = (bher_norm - 1.0 + M_LN2)/M_LN10;
  }
  else if (bher < bh_eta_crit) {
    bher_norm = 0.5;
    bher = (ledd + bh_eta_crit)/2.0;
  }
  bher-=bher_char;
  //if (bher < BHER_MIN || bher > BHER_EFF_MAX) return 0;
  double retval = bher_norm*bh_prob_norm / (exp10(bher*alpha) + exp10(bher*delta));
  if (!isfinite(retval)) {
    //fprintf(stderr, "%e %e %e %e (infinite)\n", ledd, bher_char, alpha, bh_prob_norm);
    return 0;
  }
  return retval;
}

// Note here that the input Eddington ratio (ledd) is the ***kinetic*** one,
// NOT the total Eddington ratio. This is because the conversion between 
// them is not monotonic.
double _prob_of_ledd_kinetic(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit) {
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
  if (!isfinite(retval)) {
    //fprintf(stderr, "%e %e %e %e (infinite)\n", ledd, bher_char, alpha, bh_prob_norm);
    return 0;
  }
  return retval;
}

double _prob_of_ledd_linear(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit) {
  double bher = ledd;
  double bher_norm = 1;
  bher-=bher_char;
  //  if (bher < BHER_MIN || bher > BHER_EFF_MAX) return 0;
  return bher_norm*bh_prob_norm / (exp10(bher*alpha) + exp10(bher*delta));
}

double calc_qpdf_at_l_m_z(double l, double m, double z) {
  int64_t i;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = steps[step].smhm.scatter;

  double lum_i = 72.5 - 2.5*(l - 7);
  double Lx = find_Lx_at_Lbol(l);
  double psi_ctn = psi(Lx, z);
  double acc_rate_norm = -(lum_i+5.26)/2.5;
  double med_bh_m = calc_bh_at_bm(bulge_mass(m, 1.0/(1.0+z)), steps[step].smhm);
  double med_edd_r = acc_rate_norm - med_bh_m;
  double tw = 0, ld=0;
  for (i=0; i<M_BINS; i++) {
    double dm = (steps[step].log_sm[i]+steps[step].smhm.mu - m)/s1;
    double w = exp(-0.5*dm*dm)*steps[step].t[i];
    // double f_CTK = steps[step].f_CTK[i];
    tw += w;
    double eta_frac = med_edd_r - steps[step].bh_eta[i];
    //printf("%f %f %f %f %e %f %f %f\n", qi->l, m, mbh, eta, nd, alpha, eta_0, eta_frac);
    if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) continue;
    double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
    int64_t bher_b = bher_f;
    bher_f -= bher_b;
    //double p1 = steps[step].bher_dist[bher_b];
    //double p2 = steps[step].bher_dist[bher_b+1]; 

    double p1 = steps[step].bher_dist_full[i*BHER_BINS+bher_b];
    double p2 = steps[step].bher_dist_full[i*BHER_BINS+bher_b+1];
    if (bher_b >= BHER_BINS-1) p2 = p1;
    double prob = p1 + bher_f*(p2-p1);
    // prob *= (1 - f_CTK);
    // mass-dependent modulation of the duty cycle
    // prob *= steps[step].smhm.bh_duty * pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
    // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
    // f_mass = f_mass < 1? f_mass : 1;
    double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    double dc = steps[step].smhm.bh_duty * f_mass;
    if (dc < 1e-4) dc = 1e-4;
    prob *= dc;
    ld += w*prob;
  }
  // move the following line into the loop above, due to the implementation of the mass-dependent modulation of the duty cycle.
  // ld *= steps[step].smhm.bh_duty;
  if (tw>0) return (ld/tw/(1+psi_ctn));
  return 1e-15;
}

double calc_qpdf_at_l_m_z_new(double lbol, double m, double z) {
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

  // Compton-thick fractions.
  double lx = find_Lx_at_Lbol(lbol);
  double psi_ctn = psi(lx, z);

  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  double qpdf = 0;


  //if (m >= 11.1) fprintf(stderr, "#SM=%.1f, sBHAR=%.1f\n#m Edd_r bh_eta eta_frac ledd_min ledd_max prob weight\n", m, sBHAR);
  for (i=0; i<MBH_BINS; i++) {
    double mbh = steps[step].bh_mass_min + (i + 0.5) * mbh_inv_bpdex;
    double dmbh = (mbh - mbh_med)/bh_scatter;
    // Note the mbh_inv_bpdex in the next line.
    double prob_mbh = exp(-0.5*dmbh*dmbh)*gauss_norm * mbh_inv_bpdex;
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
      // prob /= (1 + psi_ctn);
      // prob *= steps[step].smhm.bh_duty * pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
      // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
      // f_mass = f_mass < 1? f_mass : 1;
    }
    qpdf += prob_mbh * prob_lbol;
    
    //if (m >= 11.1) fprintf(stderr, "%.1f %.6f %.6f %.6f %.3f %.3f %.6e %.6e\n", 
      //                    M_MIN + (i + 0.5) * INV_BPDEX, med_edd_r, steps[step].bh_eta[i], 
        //                  med_edd_r - steps[step].bh_eta[i], steps[step].ledd_min[i],
          //                steps[step].ledd_max[i], prob, w);
  }
  // // move the following line into the loop above, due to the implementation of the mass-dependent modulation of the duty cycle.
  // // ld *= steps[step].smhm.bh_duty;
  // //if (m >= 11.1) fprintf(stderr, "tw=%.6e\n", tw);
  // //if (m >= 11.1) fprintf(stderr, "prob=%.6e\n\n",ld/tw);
  // if (tw>0) return (ld/tw);
  // return 1e-15;
  return qpdf / (1 + psi_ctn);
}

double calc_qpdf_at_sBHAR_m_z(double sBHAR, double m, double z) {
  int64_t i;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = steps[step].smhm.scatter;
  double tw = 0, ld=0;
  //if (m >= 11.1) fprintf(stderr, "#SM=%.1f, sBHAR=%.1f\n#m Edd_r bh_eta eta_frac ledd_min ledd_max prob weight\n", m, sBHAR);
  for (i=0; i<M_BINS; i++) {
        double dm = (steps[step].log_sm[i]+steps[step].smhm.mu - m)/s1;
    double w = exp(-0.5*dm*dm)*steps[step].t[i];
    // double Lbol = log10(2.6e35) + steps[step].log_sm[i]+steps[step].smhm.mu + sBHAR;
    double Lbol = log10(2.6e35) + m + sBHAR;
    double lum_i = 72.5 - 2.5*(Lbol - 7);
    double Lx = find_Lx_at_Lbol(Lbol);
    double acc_rate_norm = -(lum_i+5.26)/2.5;
    double med_bh_m = calc_bh_at_bm(bulge_mass(m, 1.0/(1.0+z)), steps[step].smhm);
    // double med_edd_r = acc_rate_norm - med_bh_m;
    double med_edd_r = acc_rate_norm - med_bh_m + steps[step].smhm.eta_mu;
    double psi_ctn = psi(Lx, z);
    // double f_CTK = steps[step].f_CTK[i];
    // fprintf(stdout, "z=%f, Mh=%f, SM=%f, Mbh=%f, bh_eta_0=%f, bh_eta=%f, Lbol=%f\n", z, M_MIN + (i + 0.5) * INV_BPDEX,
    //              steps[step].log_sm[i]+steps[step].smhm.mu, med_bh_m, steps[step].bh_eta[i], med_edd_r, Lbol);
    tw += w;
    double eta_frac = med_edd_r - steps[step].bh_eta[i];
    //printf("%f %f %f %f %e %f %f %f\n", qi->l, m, mbh, eta, nd, alpha, eta_0, eta_frac);
    double prob;
    if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) 
      prob = 0;
    else
    {
      double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t bher_b = bher_f;
      bher_f -= bher_b;
      /*double p1 = steps[step].bher_dist[bher_b];
        double p2 = steps[step].bher_dist[bher_b+1]; */

      double p1 = steps[step].bher_dist_full[i*BHER_BINS+bher_b];
      double p2 = steps[step].bher_dist_full[i*BHER_BINS+bher_b+1];
      if (bher_b >= BHER_BINS-1) p2 = p1;
      prob = p1 + bher_f*(p2-p1);
      prob /= (1 + psi_ctn);
      // prob *= steps[step].smhm.bh_duty * pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
      // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
      // f_mass = f_mass < 1? f_mass : 1;
    }
    
    double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    double dc = steps[step].smhm.bh_duty * f_mass;
    if (dc < 1e-4) dc = 1e-4;
    prob *= dc;
    ld += w*prob;

    //if (m >= 11.1) fprintf(stderr, "%.1f %.6f %.6f %.6f %.3f %.3f %.6e %.6e\n", 
      //                    M_MIN + (i + 0.5) * INV_BPDEX, med_edd_r, steps[step].bh_eta[i], 
        //                  med_edd_r - steps[step].bh_eta[i], steps[step].ledd_min[i],
          //                steps[step].ledd_max[i], prob, w);
  }
  // move the following line into the loop above, due to the implementation of the mass-dependent modulation of the duty cycle.
  // ld *= steps[step].smhm.bh_duty;
  //if (m >= 11.1) fprintf(stderr, "tw=%.6e\n", tw);
  //if (m >= 11.1) fprintf(stderr, "prob=%.6e\n\n",ld/tw);
  if (tw>0) return (ld/tw);
  return 1e-15;
}


double calc_qpdf_at_sBHAR_m_z_new(double sBHAR, double m, double z) {
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
  double lbol = log10(2.6e35) + m + sBHAR;
  // double lum_i = 72.5 - 2.5*(Lbol - 7);
  // Compton-thick fractions.
  double lx = find_Lx_at_Lbol(lbol);
  double psi_ctn = psi(lx, z);

  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;
  double eta_mu = steps[step].smhm.eta_mu;
  double qpdf = 0;


  //if (m >= 11.1) fprintf(stderr, "#SM=%.1f, sBHAR=%.1f\n#m Edd_r bh_eta eta_frac ledd_min ledd_max prob weight\n", m, sBHAR);
  for (i=0; i<MBH_BINS; i++) {
    double mbh = steps[step].bh_mass_min + (i + 0.5) * mbh_inv_bpdex;
    double dmbh = (mbh - mbh_med)/bh_scatter;
    // Note the mbh_inv_bpdex in the next line.
    double prob_mbh = exp(-0.5*dmbh*dmbh)*gauss_norm * mbh_inv_bpdex;
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
      // prob /= (1 + psi_ctn);
      // prob *= steps[step].smhm.bh_duty * pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
      // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
      // f_mass = f_mass < 1? f_mass : 1;
    }
    qpdf += prob_mbh * prob_lbol;
    
    //if (m >= 11.1) fprintf(stderr, "%.1f %.6f %.6f %.6f %.3f %.3f %.6e %.6e\n", 
      //                    M_MIN + (i + 0.5) * INV_BPDEX, med_edd_r, steps[step].bh_eta[i], 
        //                  med_edd_r - steps[step].bh_eta[i], steps[step].ledd_min[i],
          //                steps[step].ledd_max[i], prob, w);
  }
  // // move the following line into the loop above, due to the implementation of the mass-dependent modulation of the duty cycle.
  // // ld *= steps[step].smhm.bh_duty;
  // //if (m >= 11.1) fprintf(stderr, "tw=%.6e\n", tw);
  // //if (m >= 11.1) fprintf(stderr, "prob=%.6e\n\n",ld/tw);
  // if (tw>0) return (ld/tw);
  // return 1e-15;
  return qpdf / (1 + psi_ctn);
}


double calc_qpdf_at_eta_m_z(double eta, double m, double z) {
  int64_t i;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  if (step >= num_outputs - 1) { step = num_outputs - 2; f = 1.0; }
  double s1 = steps[step].smhm.scatter;

  //double lum_i = 72.5 - 2.5*(l - 7);
  //double acc_rate_norm = -(lum_i+5.26)/2.5;
  //double med_bh_m = calc_bh_at_bm(bulge_mass(m, 1.0/(1.0+z)), steps[step].smhm);
  //double med_edd_r = acc_rate_norm - med_bh_m;
  double tw = 0, ld=0;

  for (i=0; i<M_BINS; i++) {
    double dm = (steps[step].log_sm[i]+steps[step].smhm.mu - m)/s1;
    // double f_CTK = steps[step].f_CTK[i];
    double w = exp(-0.5*dm*dm)*steps[step].t[i];
    tw += w;
    //double eta_frac = med_edd_r - steps[step].bh_eta[i];
    double eta_frac = eta - steps[step].bh_eta[i];
    //printf("%f %f %f %f %e %f %f %f\n", qi->l, m, mbh, eta, nd, alpha, eta_0, eta_frac);
    if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) continue;
    double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
    int64_t bher_b = bher_f;
    bher_f -= bher_b;
    /*double p1 = steps[step].bher_dist[bher_b];
      double p2 = steps[step].bher_dist[bher_b+1]; */

    double p1 = steps[step].bher_dist_full[i*BHER_BINS+bher_b];
    double p2 = steps[step].bher_dist_full[i*BHER_BINS+bher_b+1];
    if (bher_b >= BHER_BINS-1) p2 = p1;
    double prob = p1 + bher_f*(p2-p1);
    // prob *= (1 - f_CTK);
    // prob *= steps[step].smhm.bh_duty * pow(log10(steps[step].bh_mass_avg[i]) / log10(steps[step].bh_mass_avg[M_BINS - 1]), steps[step].smhm.dc_beta);
    // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[step].smhm.dc_beta);
    // f_mass = f_mass < 1? f_mass : 1;
    double f_mass = exp((log10(steps[step].bh_mass_avg[i]) - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
  f_mass = f_mass / (1 + f_mass);
    double dc = steps[step].smhm.bh_duty * f_mass;
    if (dc < 1e-4) dc = 1e-4;
    prob *= dc;
    ld += w*prob;
  }
  // move the following line into the loop above, due to the implementation of the mass-dependent modulation of the duty cycle.
  // ld *= steps[step].smhm.bh_duty;
  if (tw>0) return (ld/tw);
  return 1e-15;
}

double calc_cosmic_sfr(double z) {
  int64_t step;
  double f;
  double sfr;
  calc_step_at_z(z, &(step), &(f));
  sfr = steps[step].observed_cosmic_sfr;
  if (f>-1) sfr += f*(steps[step+1].observed_cosmic_sfr-sfr);
  if (sfr < 1e-10) sfr = 1e-10;
  return sfr;
}

double calc_cosmic_bhar(double z) {
  int64_t step;
  double f;
  double bhar;
  calc_step_at_z(z, &(step), &(f));
  bhar = steps[step].observed_cosmic_bhar;
  if (f>-1) bhar += f*(steps[step+1].observed_cosmic_bhar-bhar);
  if (bhar < 1e-10) bhar = 1e-10;
  // printf("cosmic_bhar=%e\n", bhar);
  return bhar;
}

double chi2_err_helper(double Vc, void *extra_data) {
  double m = *((double *)extra_data);
  double z = comoving_volume_to_redshift(Vc);
  double smf = evaluate_from_step(z, m, 1);
  if (smf < 1e-15) smf = 1e-15;
  return smf;
}

double chi2_err_helper_qf(double Vc, void *extra_data) {
  double m = *((double *)extra_data);
  double z = comoving_volume_to_redshift(Vc);
  double qf = evaluate_from_step(z, m, 2);
  // if (smqff < 1e-15) smf = 1e-15;
  return qf;
}

double chi2_err_helper_uv(double Vc, void *extra_data) {
  double uv = *((double *)extra_data);
  double z = comoving_volume_to_redshift(Vc);
  double uvlf = evaluate_from_step_uv(z, uv);
  if (uvlf < 1e-15) uvlf = 1e-15;
  return uvlf;
}

//inline 
double model_chi2(struct obs_smf *m_smf,
			struct obs_smf *r_smf) {
  int i;
  double chi2 = 0;
  double err, real_err;
  struct real_smf *real_smf = r_smf->real_smf;
  struct real_smf *model_smf = m_smf->real_smf;
  for (i=0; i<r_smf->real_smf_points; i++) {
    double fit_error = (r_smf->linear_errors) ? fabs(real_smf[i].val*(exp10fc(FIT_TOLERANCE)-1.0)) : fabs(FIT_TOLERANCE);
    err = (model_smf[i].val-real_smf[i].val);
    if (fabs(err) <= fit_error) continue;
    if (err > 0) err -= fit_error;
    else err += fit_error;

    if (err > real_smf[i].err_h) err/=real_smf[i].err_h;
    else if (-err > real_smf[i].err_l) err/=real_smf[i].err_l;
    else {
      real_err = real_smf[i].err_l +
	(err+real_smf[i].err_l)/(real_smf[i].err_h+real_smf[i].err_l) *
	(real_smf[i].err_h - real_smf[i].err_l);
      err /= real_err;
    }
    if (!isfinite(err)) {
     fprintf(stderr, "Infinite error at z=%f, m=%f, type=%d, model val=%f, real val=%f\n", 
	      0.5*(r_smf->z_low + r_smf->z_high),
	      real_smf[i].mass, r_smf->type, model_smf[i].val, real_smf[i].val);
      // chi2 = 1e30;
      // return chi2;
    }
    chi2 += err*err;
  }
  return chi2;
}

//inline 
double model_chi2_point(struct obs_smf_point *m_smf,
			      struct obs_smf_point *r_smf) {
  double chi2 = 0;
  double err, real_err;
  double fit_error = (r_smf->linear_errors) ? fabs(r_smf->val*(pow(10,FIT_TOLERANCE)-1.0)) : fabs(FIT_TOLERANCE);
  err = (m_smf->val-r_smf->val);
  if (fabs(err) <= fit_error) 
  {
   //if (r_smf->type == SMF_TYPE)   
     // fprintf(stdout, "%f %f %f %f %f %f\n",
     //r_smf->z_low, r_smf->z_high, r_smf->mass, r_smf->val, m_smf->val, 0.0);	
      return 0;
  }
  if (err > 0) err -= fit_error;
  else err += fit_error;

  if (err > r_smf->err_h) err/=r_smf->err_h;
  else if (-err > r_smf->err_l) err/=r_smf->err_l;
  else {
    real_err = r_smf->err_l +
      (err+r_smf->err_l)/(r_smf->err_h+r_smf->err_l) *
      (r_smf->err_h - r_smf->err_l);
    err /= real_err;
  }
  //if (r_smf->type == SMF_TYPE) 
    // fprintf(stderr, "Type %d at z_low=%f, z_high=%f, mass=%f, data=%f, model=%f, chi2=%f\n", r_smf->type,
     //fprintf(stdout, "%f %f %f %f %f %f\n",
     //r_smf->z_low, r_smf->z_high, r_smf->mass, r_smf->val, m_smf->val, err * err);
  if (!isfinite(err)) {
    fprintf(stderr, "Infinite error at z=%f, m=%f, type=%d, model val=%f, real val=%f\n", 
	    0.5*(r_smf->z_low + r_smf->z_high),
	    r_smf->mass, r_smf->type, m_smf->val, r_smf->val);
    // chi2 = 1e30;
    //   return chi2;
  }
  chi2 += err*err;
  return chi2;
}

double calc_single_chi2_err(struct obs_smf *cur_smf) {
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

  for (i=0; i<real_smf_points; i++) {
    m = real_smf[i].mass;
    if (cur_smf->type == SMF_TYPE) {
      if (cur_smf->z_low != cur_smf->z_high && !no_z_scaling) {
	if (cur_smf->z_low < 0.2) {
	   //epsilon *= 0.01;
	   smf_val = adaptiveGauss(chi2_err_helper, &m, v_low, v_high,
	 			  PHI_INTEGRAL_PRECISION*0.01, omp_get_thread_num()+20);
	 }
	 else {
	  epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
	  smf_val = adaptiveSimpsons(chi2_err_helper, &m,
				     v_low, v_high, epsilon, 10);
	 }
	smf_val /= weight;
      }
      else {
	smf_val = chi2_err_helper(0.5*(v_low+v_high), &m);
      }
    }

    else if (cur_smf->type == QF_TYPE) 
    {
        // if (cur_smf->z_low != cur_smf->z_high && !no_z_scaling) 
        // {
        //   if (cur_smf->z_low < 0.2) 
        //   {
        //     //epsilon *= 0.01;
        //     smf_val = adaptiveGauss(chi2_err_helper_qf, &m, v_low, v_high,
        //           PHI_INTEGRAL_PRECISION*0.01, omp_get_thread_num()+20);
        //   }
        //   else 
        //   {
        //     epsilon = chi2_err_helper_qf((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
        //     smf_val = adaptiveSimpsons(chi2_err_helper_qf, &m,
        //              v_low, v_high, epsilon, 10);
        //   }
        //   smf_val /= weight;
        // }
        // else 
        // {
      smf_val = chi2_err_helper_qf(0.5*(v_low+v_high), &m);
        // }
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
      smf_val = calc_cosmic_sfr(real_smf[i].mass); //Actually redshift for CSFR
    } 
    // else if (cur_smf->type == QLF_TYPE) 
    // {
    //   smf_val = calc_quasar_lf(m, cur_smf->z_low);
    // }
    // else if (cur_smf->type == QLF_CTN_TYPE) 
    // {
    //   smf_val = calc_quasar_lf_ctn(m, cur_smf->z_low);
    // }
    // else if (cur_smf->type == QPDF_TYPE) 
    // {
    //   smf_val = calc_qpdf_at_l_m_z(cur_smf->extra, cur_smf->mass, median_z);
    // }
    // else if (cur_smf->type == QPDF_ETA_TYPE)
    // {
    //   smf_val = calc_qpdf_at_sBHAR_m_z(cur_smf->extra, cur_smf->mass, median_z);
    // }
    // else //if (cur_smf->type == ABHMF_TYPE) 
    // {
    //   smf_val = calc_active_bhmf(cur_smf->mass, median_z);
    // }

    else if (cur_smf->type == BHMF_TYPEI_TYPE)
    {
      smf_val = calc_bhmf_typeI(m, cur_smf->z_low);
    }



    model_smf.real_smf[i].mass = m;
    model_smf.real_smf[i].val = (cur_smf->linear_errors)
      ? smf_val : log10(smf_val);
    //printf("Z=%f-%f; M: %f; Type: %d; Val: %e\n", cur_smf->z_low, cur_smf->z_high, m, cur_smf->type, smf_val);
  }
  chi2 = model_chi2(&model_smf, cur_smf);
  return chi2;
}

double calc_single_chi2_err_point(struct obs_smf_point *cur_smf) {
  double smf_val, epsilon;
  struct obs_smf_point model_smf;
  double v_high = comoving_volume(cur_smf->z_high);
  double v_low = comoving_volume(cur_smf->z_low);
  double weight = v_high - v_low;
  double m;
  double median_z = cbrt(0.5*(cur_smf->z_high*cur_smf->z_high*cur_smf->z_high + cur_smf->z_low*cur_smf->z_low*cur_smf->z_low));


  m = cur_smf->mass;
  if (cur_smf->type == SMF_TYPE) {
    if (cur_smf->z_low != cur_smf->z_high) {
 //      if (cur_smf->z_low < 0.2) {
	// //epsilon *= 0.01;
	// smf_val = adaptiveGauss(chi2_err_helper, &m, v_low, v_high,
	// 			PHI_INTEGRAL_PRECISION*0.01, omp_get_thread_num()+10);
 //      }
 //      else {
	epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
	smf_val = adaptiveSimpsons(chi2_err_helper, &m,
				   v_low, v_high, epsilon, 10);
 //      }
      smf_val /= weight;
      // smf_val = chi2_err_helper((v_high+v_low)/2.0, &m);
    }
    else {
      smf_val = chi2_err_helper(v_low, &m);
    }
  }
  else if (cur_smf->type == SSFR_TYPE) {
    smf_val = calc_ssfr(m, cur_smf->z_low);

  }
  // else { //if (cur_smf->type == CSFR_TYPE) {
  //   smf_val = calc_cosmic_sfr(cur_smf->mass); //Actually redshift for CSFR
  // }
  
  else if (cur_smf->type == CSFR_TYPE) {
    smf_val = calc_cosmic_sfr(cur_smf->mass); //Actually redshift for CSFR
  } 
  else if (cur_smf->type == QLF_TYPE) 
  {
    smf_val = calc_quasar_lf_new(m, cur_smf->z_low);
  }
  else if (cur_smf->type == QLF_CTN_TYPE) 
  {
    smf_val = calc_quasar_lf_ctn(m, cur_smf->z_low);
  }
  else if (cur_smf->type == QPDF_TYPE) 
  {
    smf_val = calc_qpdf_at_l_m_z_new(cur_smf->extra, cur_smf->mass, median_z);
    //smf_val = calc_qpdf_at_l_m_z(cur_smf->extra, cur_smf->mass, median_z);
  }
  else if (cur_smf->type == QPDF_ETA_TYPE)
  {
    smf_val = calc_qpdf_at_sBHAR_m_z_new(cur_smf->extra, cur_smf->mass, median_z);
  }
  else if (cur_smf->type == ABHMF_TYPE) 
  {
    smf_val = calc_active_bhmf(cur_smf->mass, median_z);
  }

  else if (cur_smf->type == QF_TYPE) {
  //   if (cur_smf->z_low != cur_smf->z_high) {
  //     if (cur_smf->z_low < 0.2) {
  // //epsilon *= 0.01;
  // smf_val = adaptiveGauss(chi2_err_helper_qf, &m, v_low, v_high,
  //       PHI_INTEGRAL_PRECISION*0.01, omp_get_thread_num()+10);
  //     }
  //     else {
  // epsilon = chi2_err_helper_qf((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
  // smf_val = adaptiveSimpsons(chi2_err_helper_qf, &m,
  //          v_low, v_high, epsilon, 10);
  //     }
  //     smf_val /= weight;
  //   }
    // else {
      smf_val = chi2_err_helper_qf(0.5*(v_low + v_high), &m);
    // }
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
  model_smf.val = (cur_smf->linear_errors)
    ? smf_val : log10(smf_val);
  //if (cur_smf->type == QLF_TYPE)
  //{
    //printf("%f %f %e %f\n", 0.5 * (cur_smf->z_low + cur_smf->z_high), model_smf.mass, model_smf.val, model_chi2_point(&model_smf, cur_smf));
  //}
  return model_chi2_point(&model_smf, cur_smf);
}


//inline 
double _mf_cache(struct mf_cache *cache, double mass) {
  int i;
  double m, r;
  if (mass >= cache->mass_max) {
    return -1000;
    //m = (mass - cache->mass_min);
    //r = (cache->mf[0] + (1.0-cache->alpha)*(m) - exp10fc(mass+cache->inv_m0));
  }
  else {
    m = (mass - cache->mass_min)*cache->inv_mass_spacing;
    if (m < 0) r=((1.0-cache->alpha)*(mass-cache->mass_min) + cache->mf[0]);
    else {
      i = m;
      r = (cache->mf[i] + (cache->mf[i+1]-cache->mf[i])*(m-i));
    }
  }
  return r;
}

double mf_cache(double scale, double mass) {
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


void load_real_smf(struct obs_smf_point **cur_smf, int64_t *num_points, char *filename)
{
  FILE *smf_data;
  char buffer[1024];
  // char *data_types[] = {"smf", "ssfr", "cosmic sfr"};
  char *data_types[] = {"smf", "ssfr", "cosmic sfr", "quasar lf", "quasar luminosity pdf", "active bhmf", "quasar eddington ratio pdf", "Compton-thin quasar lf", "qf", "uvlf"};
  char *error_types[] = {"log", "linear"};
  struct obs_smf_point next_smf, conv_smf;
  int64_t n, loaded=0;
  smf_data = fopen(filename, "r");
  if (!smf_data) {
    printf("Couldn't open file %s for reading!\n", filename);
    exit(1);
  }

  next_smf.type = SMF_TYPE;
  next_smf.linear_errors = 0;
  next_smf.z_low = next_smf.z_high = 0;
  while (fgets(buffer, 1023, smf_data)) {
    n = sscanf(buffer, "%lf %lf %lf %lf", &(next_smf.mass),
	       &(next_smf.val), &(next_smf.err_h), &(next_smf.err_l));
    if (buffer[0] == '#') {
      if (!strncmp(buffer, "#zlow: ", 7)) next_smf.z_low = atof(&(buffer[7]));
      if (!strncmp(buffer, "#zhigh: ", 8)) next_smf.z_high = atof(&(buffer[8]));
      if (!strncmp(buffer, "#errors: ", 9)) {
	if (!strncasecmp(&(buffer[9]), "Linear", 6)) next_smf.linear_errors = 1;
	else  next_smf.linear_errors = 0;
      }
      if (!strncmp(buffer, "#type: ", 7)) {
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
    if (n!=4 || !next_smf.err_l || !next_smf.err_h) {
      if (strlen(buffer) > 1)
	fprintf(stderr, "Invalid line (skipping): %s\n", buffer);
      continue;
    }

    if (next_smf.type == QPDF_TYPE || next_smf.type == QPDF_ETA_TYPE) {
      n = sscanf(buffer, "%lf %lf %lf %lf %lf", &(next_smf.mass), &(next_smf.extra),
         &(next_smf.val), &(next_smf.err_h), &(next_smf.err_l));

      if (n!=5 || !next_smf.err_l || !next_smf.err_h) {
  if (strlen(buffer) > 1)
    fprintf(stderr, "Invalid line (skipping): %s\n", buffer);
  continue;
      }
    }
    
    //Add in systematic errors from calculations.
    double syst_err = CALCULATION_TOLERANCE;
    conv_smf = next_smf;
    if (conv_smf.linear_errors) {
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
    if (!((*num_points)%100)) {
      *cur_smf = check_realloc(*cur_smf, sizeof(struct obs_smf_point)*
			       ((*num_points)+100), "Allocating SMFs.");
    }

    cur_smf[0][*num_points] = conv_smf;
    *num_points = (*num_points) + 1;
    loaded++;
  }
  fclose(smf_data);
  if (loaded) {
    if (next_smf.z_high > z_max) z_max = next_smf.z_high;
    if (next_smf.z_low < z_min) z_min = next_smf.z_low;
  }

  fprintf(stderr, "#Loaded %"PRId64" points from %s (type: %s; errors: %s)\n", loaded, filename, data_types[next_smf.type], error_types[next_smf.linear_errors]); 
}

void setup_psf(int scatter) {
  init_gauss_cache();
  use_obs_psf = scatter;
}

float readfloat(FILE *input) {
  char buffer[50] = {0};
  float result;
  fgets(buffer, 50, input);
  sscanf(buffer, "%f", &result);
  return result;
}

void load_mf_cache(char *filename) {
  FILE *input;
  int i;
  float omega_m, omega_l, h0;
  struct mf_cache *mf;
  if (!(input = fopen(filename, "r"))) {
    printf("Couldn't open MF cache %s!\n", filename);
    exit(6);
  }

  all_mf_caches = (struct z_mf_cache *)malloc(sizeof(struct z_mf_cache));

  omega_m = all_mf_caches->omega_m = readfloat(input);
  omega_l = all_mf_caches->omega_l = readfloat(input);
  h0 = all_mf_caches->h0 = readfloat(input);
  init_cosmology(omega_m, omega_l, h0);
  init_time_table(omega_m, h0);

  all_mf_caches->scale_min = readfloat(input);
  all_mf_caches->scale_max = readfloat(input);
  all_mf_caches->scales = readfloat(input);
  all_mf_caches->avg_inv_scale_spacing = (float)(all_mf_caches->scales-1) / 
    (all_mf_caches->scale_max - all_mf_caches->scale_min);

  if (! (all_mf_caches->scale_min > 0 && all_mf_caches->scale_min < 1.05 &&
	 all_mf_caches->scale_max > all_mf_caches->scale_min &&
	 all_mf_caches->scale_max < 1.05)) {
    printf("Invalid scales in MF cache (%f - %f), expecting 0<a<1.\n",
	   all_mf_caches->scale_min, all_mf_caches->scale_max);
    exit(7);
  }

  all_mf_caches->caches = (struct mf_cache *)
    malloc(sizeof(struct mf_cache)*all_mf_caches->scales);
  for (i=0; i<all_mf_caches->scales; i++) {
    mf = &(all_mf_caches->caches[i]);
    mf->mass_min = readfloat(input);
    mf->mass_max = readfloat(input);
    mf->masses = readfloat(input);
    mf->scale = readfloat(input);
    mf->alpha = fabs(readfloat(input));
    if (! (mf->scale >= all_mf_caches->scale_min &&
	   mf->scale <= all_mf_caches->scale_max)) {
      printf("Invalid scale in MF cache index %d: %f\n", i, mf->scale);
      exit(7);
    }
    if (! (mf->alpha > 1 && mf->alpha < 3)) {
      printf("Invalid |alpha| in MF cache: %f (Expected 1 - 3)\n", mf->alpha);
      exit(7);
    }

    mf->mf = (float *)malloc(sizeof(float)*mf->masses);
    fread(mf->mf, sizeof(float), mf->masses, input);
    mf->inv_mass_spacing = (float)(mf->masses-1)/(mf->mass_max - mf->mass_min);
    mf->inv_m0 = log10(((mf->mass_max - mf->mass_min)*(1.0-mf->alpha)
			+ mf->mf[0] - mf->mf[mf->masses-1])/exp10(mf->mass_max));
  }
  for (i=0; i<all_mf_caches->scales-1; i++)
    all_mf_caches->caches[i].inv_scale_spacing = 
      1.0 / (all_mf_caches->caches[i+1].scale - all_mf_caches->caches[i].scale);
  fclose(input);

  gen_exp10cache(); //Required for using MF cache
}
