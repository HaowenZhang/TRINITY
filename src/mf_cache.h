#ifndef MF_CACHE
#define MF_CACHE

struct mf_cache {
  float mass_min;
  float mass_max;
  float inv_scale_spacing;
  float inv_mass_spacing;
  int masses;
  float scale;
  float alpha;
  float inv_m0;
  float *mf;
  float *errors;
};

struct z_mf_cache {
  float scale_min;
  float scale_max;
  float h0, omega_m, omega_l;
  float avg_inv_scale_spacing;
  int scales;
  struct mf_cache *caches;
};

#endif /* MF_CACHE */
