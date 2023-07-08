#ifndef _INTHASH_H_
#define _INTHASH_H_
#include <inttypes.h>

#define IH_INVALID INT64_MAX
#define IH_DELETED (INT64_MAX-1)
#define IH_SKIP (INT64_MAX-2)

struct intbucket {
  int64_t key;
  void *data;
};

struct inthash {
  uint64_t hashwidth, elems, num_buckets, hashnum;
  struct intbucket *buckets;
};


struct inthash *new_inthash(void);
void *ih_keylist(struct inthash *ih);
void *ih_getval(struct inthash *ih, int64_t key);
void ih_setval(struct inthash *ih, int64_t key, void *data);
void ih_delval(struct inthash *ih, int64_t key);
void free_inthash(struct inthash *ih);
void ih_prealloc(struct inthash *ih, int64_t size);

void free_inthash2(struct inthash *ih);
void ih_setval2(struct inthash *ih, int64_t key1, int64_t key2, void *data);
void *ih_getval2(struct inthash *ih, int64_t key1, int64_t key2);
void ih_setint64(struct inthash *ih, int64_t key, int64_t value);
int64_t ih_getint64(struct inthash *ih, int64_t key);

#endif /* _INTHASH_H_ */

