#ifndef MT_RAND_H
#define MT_RAND_H


#ifdef __cplusplus
extern "C" {
#endif

void init_genrand64(unsigned long long seed);
double genrand64_real2(void);
float normal_random(float mean, float stddev);
unsigned long long genrand64_int64(void);

//R250 Compatibility routines
#define r250_init init_genrand64
#define r250 genrand64_int64
#define dr250 genrand64_real2

#ifdef __cplusplus
}
#endif


#endif /* MT_RAND_H */
