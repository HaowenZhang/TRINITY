#ifndef PSF_CACHE
#define PSF_CACHE

#define PSF_MASS_MIN -4
#define PSF_MASS_MAX 4
#define PSF_CACHE_SIZE 1024
#define PSF_MASS_STEP ((PSF_MASS_MAX-PSF_MASS_MIN)/((double)PSF_CACHE_SIZE))
#define PSF_INV_MASS_STEP (1.0/PSF_MASS_STEP)
#define PSF_SCATTER_MIN 0.0
#define PSF_SCATTER_MAX 0.5
#define PSF_SCATTER_STEPS 100
#define PSF_SCATTER_STEP ((PSF_SCATTER_MAX-PSF_SCATTER_MIN)/ \
			  ((double)PSF_SCATTER_STEPS))
#define PSF_SCATTER_EFF_MAX (PSF_SCATTER_MAX-1.05*PSF_SCATTER_STEP)
#define PSF_INV_SCATTER_STEP (1.0/PSF_SCATTER_STEP)
#define PSF_OBS_MIN 0.0
#define PSF_OBS_MAX 0.6
#define PSF_OBS_STEPS 120
#define PSF_OBS_STEP ((PSF_OBS_MAX-PSF_OBS_MIN)/((double)PSF_OBS_STEPS))
#define PSF_INV_OBS_STEP (1.0/PSF_OBS_STEP)
#define PSF_OBS_EFF_MAX (PSF_OBS_MAX-2.1*PSF_OBS_STEP)

struct psf_cache {
  float scatter;
  float obs;
  float min_mass, max_mass;
  float cache[PSF_CACHE_SIZE];
};

#endif /* PSF_CACHE */
