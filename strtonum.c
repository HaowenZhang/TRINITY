#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#ifdef CONV_64BIT

#define FLOAT_TYPE double
#define FLOAT_LONGER_TYPE long double 
#define FLOAT_MAX DBL_MAX
#define FLOAT_EXPS DBL_MAX_10_EXP
#define FLOAT_CONV_NAME strtodouble
#define ITYPE int64_t
#define IMAX INT64_MAX
#define IMIN INT64_MIN
#define ICONV_NAME strtoint64

#else 

#define FLOAT_TYPE float
#define FLOAT_MAX FLT_MAX
#define FLOAT_EXPS FLT_MAX_10_EXP
#define FLOAT_CONV_NAME strtofloat
#define FLOAT_LONGER_TYPE double
#define ITYPE int32_t
#define IMAX INT32_MAX
#define IMIN INT32_MIN
#define ICONV_NAME strtoint32

#endif /* CONV_64BIT */


static inline FLOAT_TYPE FLOAT_CONV_NAME (const char *cur_pos, char **end_pos)
{
  FLOAT_LONGER_TYPE sign = 1;
  FLOAT_LONGER_TYPE val = 0;
  FLOAT_LONGER_TYPE savedval = 0;
  const char *periodpos=0, *periodendpos=0;
  const char *pos = cur_pos;
  const char zero = '0';
  const char zerom1 = '0'-1;
  const char ninep1 = '9'+1;
  const FLOAT_LONGER_TYPE toobig = FLOAT_MAX / 20.0;
  int seenexp = 0;
  int seenexpsign = 0;
  char c;
  static int generated_exps = 0;
  static FLOAT_LONGER_TYPE exps[FLOAT_EXPS*4 + 1];
  int i, e = 0;

  if (!generated_exps) {
    generated_exps = 1;
    for (i=0; i<(FLOAT_EXPS*4+1); i++) {
      exps[i] = powl(10, i-(FLOAT_EXPS*2));
    }
  }

  if (*pos == '+') pos++;
  else if (*pos == '-') { sign = -1; pos++; }
  for (; *pos; pos++) {
    c = *pos;
    if ((c > zerom1) && (c < ninep1)) {
      if (val < toobig)
	val = val*10.0 + (c-zero);
      else
	e++;
    }
    else if (c == '.') {
      if (seenexp) break;
      periodpos = pos;
    }
    else if (c == 'e' || c == 'E') {
      if (seenexp) break;
      savedval = val;
      seenexp = 1;
      val = 0;
      periodendpos = pos;
    }
    else if (c == '+') {
      if (!seenexp || seenexpsign || ((pos - periodendpos)>1)) break;
      seenexpsign = 1;
    }
    else if (c == '-') {
      if (!seenexp || seenexpsign  || ((pos - periodendpos)>1)) break;
      seenexpsign = -1;
    }
    else break;
  }
  if (end_pos) *end_pos = (char *)pos;

  if (!seenexp) {
    savedval = val;
  } else {
    if (seenexpsign) val = copysign(val, seenexpsign);
    e += val;
  }

  if (periodpos) {
    if (!periodendpos) periodendpos = pos;
    e -= periodendpos - periodpos - 1;  
  }

  if (savedval==0) return copysign(0, sign);

  if (e > (FLOAT_EXPS*2)) {
    if (sign > 0) return FLOAT_MAX;
    else return -FLOAT_MAX;
  }
  else if (e < -(FLOAT_EXPS*2))
    return copysign(0,sign);
  
  return (copysign(savedval*exps[e+(FLOAT_EXPS*2)], sign));
}

static inline ITYPE ICONV_NAME (const char *cur_pos, char **end_pos)
{
  int sign = 1;
  ITYPE val = 0;
  const char *pos = cur_pos;
  const char zero = '0';
  const char zerom1 = '0'-1;
  const char ninep1 = '9'+1;
  const ITYPE toobig = (IMAX / 10);
  const char *exppos = 0;
  char c;
  int toobig_flag = 0;
  int seen_period = 0, seen_exp = 0;

  if (*pos == '+') pos++;
  else if (*pos == '-') { sign = -1; pos++; }
  for (; *pos; pos++) {
    c = *pos;
    if ((c > zerom1) && (c < ninep1)) {
      if (val < toobig)
	val = val*10 + (c-zero);
      else
	toobig_flag = 1;
    }
    else if (c == '.') {
      seen_period = 1;
      pos++;
      break;
    }
    else if (c == 'e' || c == 'E') {
      seen_exp = 1;
      exppos = pos;
      pos++;
      break;
    }
    else break;
  }

  if (seen_period || seen_exp) {
    for (; *pos; pos++) {
      c = *pos;
      if ((c > zerom1) && (c < ninep1)) continue;
      if (c == 'e' || c == 'E') {
	if (seen_exp) break;
	seen_exp = 1;
	exppos = pos;
	continue;
      }
      if (c == '+' || c == '-') {
	if (!seen_exp || (pos - exppos)>1) break;
	continue;
      }
      break;
    }
  }

  if (end_pos) *end_pos = (char *)pos;

  if (toobig_flag)
    return ((sign>0) ? IMAX : IMIN);

  return sign*val;
}

#undef FLOAT_TYPE
#undef FLOAT_LONGER_TYPE
#undef FLOAT_MAX
#undef FLOAT_EXPS
#undef FLOAT_CONV_NAME
#undef ITYPE
#undef IMAX
#undef IMIN
#undef ICONV_NAME

