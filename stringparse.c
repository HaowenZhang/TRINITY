#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "stringparse.h"

#include "strtonum.c"
#define CONV_64BIT
#include "strtonum.c"

static inline int strtowhatever(char *cur_pos, char **end_pos, void *data, enum parsetype type) {
  float valf;
  double vald;
  int32_t vali32;
  int64_t vali64;
  char *new_pos = cur_pos;
  char *data_pos = (char *)data;
  switch (type) {
  case PARSE_FLOAT32:
    valf = strtofloat(cur_pos, end_pos);
    if (cur_pos == *end_pos) return 0;
    *((float *)data) = valf;
    return 1;

  case PARSE_INT32:
    vali32 = strtoint32(cur_pos, end_pos);
    if (cur_pos == *end_pos) return 0;
    *((int32_t *)data) = vali32;
    return 1;

  case PARSE_FLOAT64:
    vald = strtodouble(cur_pos, end_pos);
    if (cur_pos == *end_pos) return 0;
    *((double *)data) = vald;
    return 1;

  case PARSE_INT64:
    vali64 = strtoint64(cur_pos, end_pos);
    if (cur_pos == *end_pos) return 0;
    *((int64_t *)data) = vali64;
    return 1;
   
  case PARSE_STRING:
    while (*new_pos && !(*new_pos==' ' || *new_pos=='\t' || *new_pos=='\n')) {
      *data_pos = *new_pos;
      data_pos++;
      new_pos++;
    }
    if (new_pos==cur_pos) return 0;
    *data_pos = 0;
    *end_pos = new_pos;
    return 1;

  case PARSE_SKIP:
    while (*new_pos && !(*new_pos==' ' || *new_pos=='\t' || *new_pos=='\n'))
      new_pos++;
    if (new_pos==cur_pos) return 0;
    *end_pos = new_pos;
    return 1;
  }
  return 0;
}

int64_t stringparse(char *buffer, void **data, enum parsetype *types, int64_t max_n) {
  int64_t num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  if (max_n < 1) return 0;

  while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
  while ((*cur_pos) && (num_entries < max_n)) {
   if (!strtowhatever(cur_pos, &end_pos, data[num_entries], types[num_entries]))
      break;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
  }
  return num_entries;
}

int64_t stringparse_format(char *buffer, struct parse_format *pf, int64_t max_n)
{
  int64_t num_entries = 0, np = 0;
  char *cur_pos = buffer, *end_pos;
  if (max_n < 1) return 0;
  for (np=max_n-1; np>0; np--) assert(pf[np].index > pf[np-1].index);

  while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
  while ((*cur_pos) && (np < max_n)) {
    if (num_entries == pf[np].index) {
      if (!strtowhatever(cur_pos, &end_pos, pf[np].data, pf[np].type)) break;
      np++;
    } else {
      if (!strtowhatever(cur_pos, &end_pos, NULL, PARSE_SKIP)) break;
    }
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
  }
  return np;
}

/*
int main(void) {
  int64_t i;
  double d=0;
  int d2;
  printf("%f\n", strtofloat("foo", NULL));
  printf("%f %f\n", strtofloat("1", NULL), 1.0);
  printf("%f\n", strtofloat("2", NULL));
  printf("%f\n", strtofloat("0.02", NULL));
  printf("%f\n", strtofloat("-0.02", NULL));
  printf("%f\n", strtofloat("-0.02e10foo", NULL));
  printf("%e\n", strtofloat("-0.02e-10", NULL));
  printf("%e\n", strtofloat("nAn", NULL));
  printf("%e\n", strtodouble("-inf", NULL));
  printf("%.10e\n", strtof("-0.123456789", NULL));
  printf("%.10e\n", strtofloat("-0.123456789", NULL));
  printf("%d\n", strtoint32("5002342", NULL));
  printf("%d\n", strtoint32("-3234.8345e10", NULL)); 
  printf("%ld\n", strtol("-3234.8345e10", NULL, 10)); 
  sscanf("-3234.2345e10", "%d", &d2);
  printf("%d\n", d2); 

  
  for (i=0; i<1e7; i++) {
    d += strtoint64("100.2341234", NULL);//, NULL, 10);
  }

  return 0;
}

*/
