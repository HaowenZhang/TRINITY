#ifndef _STRINGPARSE_H_
#define _STRINGPARSE_H_
#include <inttypes.h>

enum parsetype {
  PARSE_FLOAT32 = 0,
  PARSE_INT32 = 1,
  PARSE_FLOAT64 = 2,
  PARSE_INT64 = 3,
  PARSE_STRING = 4,
  PARSE_SKIP = 5,
};

#define short_parsetype parsetype

#define SHORT_PARSETYPE \
  const enum parsetype F = PARSE_FLOAT32; (void)F;	\
  const enum parsetype D = PARSE_INT32;	(void)D;	\
  const enum parsetype F64 = PARSE_FLOAT64; (void)F64;	\
  const enum parsetype LF = PARSE_FLOAT64; (void)LF;	\
  const enum parsetype LD = PARSE_INT64; (void)LD;	\
  const enum parsetype D64 = PARSE_INT64; (void)D64;	\
  const enum parsetype S = PARSE_STRING; (void)S;	\
  const enum parsetype K = PARSE_SKIP;	(void)K;	

struct parse_format {
  int64_t index;
  enum parsetype type;
  void *data;
};

int64_t stringparse(char *buffer, void **data, enum parsetype *types, int64_t max_n);
int64_t stringparse_format(char *buffer, struct parse_format *pf, int64_t max_n);

#endif
