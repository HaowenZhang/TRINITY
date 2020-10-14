#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include "check_syscalls.h"
#include "inthash.h"

#define MAX_LOAD_FACTOR 0.7

struct inthash *new_inthash(void) {
  struct inthash *ih = check_realloc(NULL, sizeof(struct inthash),
				       "Allocating inthash.");
  int64_t i;
  memset(ih, 0, sizeof(struct inthash));
  ih->hashnum = rand() + 
    (((uint64_t)rand())<<(uint64_t)32) + (uint64_t)(rand());
  ih->hashwidth = 8;
  ih->num_buckets = (uint64_t)1 << ih->hashwidth;
  ih->buckets = check_realloc(NULL, sizeof(struct intbucket)*ih->num_buckets,
			      "Allocating hash buckets.");
  memset(ih->buckets, 0, sizeof(struct intbucket)*ih->num_buckets);
  for (i=0; i<ih->num_buckets; i++) ih->buckets[i].key = IH_INVALID;
  if (!(ih->hashnum & 1)) ih->hashnum++;
  return ih;
}

void *ih_keylist(struct inthash *ih) {
  int64_t i, j=0;
  int64_t *kl = check_realloc(NULL, sizeof(int64_t)*ih->elems, "Allocating key list.");
  struct intbucket *ib = NULL;
  for (i=0; i<ih->num_buckets; i++) {
    ib = ih->buckets + i;
    if (ib->key>IH_SKIP) continue;
    kl[j]=ib->key;
    j++;
  }
  return kl;
}


inline uint64_t _ih_hash_function(struct inthash *ih, uint64_t key) {
  return ((key*ih->hashnum)>>(64 - ih->hashwidth));
}

void _ih_add_more_buckets(struct inthash *ih, int64_t add_factor) {
  int64_t i;
  int64_t old_num_buckets = ih->num_buckets;
  struct intbucket *new_buckets, *ib, *newb;

  ih->num_buckets <<= add_factor;
  int64_t new_alloc_size = sizeof(struct intbucket)*(ih->num_buckets);

  new_buckets = check_realloc(NULL, new_alloc_size, "Allocating new buckets.");
  memset(new_buckets, 0, new_alloc_size);
  for (i=0; i<ih->num_buckets; i++) new_buckets[i].key = IH_INVALID;

  ih->hashwidth+=add_factor;
  for (i=0; i<old_num_buckets; i++) {
    ib = ih->buckets + i;
    if (ib->key>IH_SKIP) continue;
    int64_t key = ib->key;
    do {
      newb = new_buckets + _ih_hash_function(ih, key);
      key = key*key + key + 3;
    } while (newb->key!=IH_INVALID);
    newb->key = ib->key;
    newb->data = ib->data;
  }

  free(ih->buckets);
  ih->buckets = new_buckets;
}


void ih_prealloc(struct inthash *ih, int64_t size) {
  if (size <= 0) return;
  int64_t numbits = ceil(log(size/MAX_LOAD_FACTOR)/log(2));
  if (numbits <= ih->hashwidth) return;
  _ih_add_more_buckets(ih, numbits-ih->hashwidth);
}

inline struct intbucket *_ih_getval(struct inthash *ih, int64_t key) {
  struct intbucket *ib = ih->buckets + _ih_hash_function(ih, key);
  int64_t key2 = key;
  while (ib->key!=IH_INVALID) {
    if (ib->key == key) return ib;
    key2 = key2*key2 + key2 + 3;
    ib = ih->buckets + _ih_hash_function(ih, key2);
  }
  return ib;
}

inline struct intbucket *_ih_getval_deleted(struct inthash *ih, int64_t key) {
  struct intbucket *ib = ih->buckets + _ih_hash_function(ih, key);
  struct intbucket *ib_del = NULL;
  int64_t key2 = key;
  while (ib->key!=IH_INVALID) {
    if (ib->key == key) return ib;
    if (ib->key == IH_DELETED) ib_del = ib;
    key2 = key2*key2 + key2 + 3;
    ib = ih->buckets + _ih_hash_function(ih, key2);
  }
  if (ib_del) return ib_del;
  return ib;
}

void *ih_getval(struct inthash *ih, int64_t key) {
  struct intbucket *ib = _ih_getval(ih, key);
  if (ib->key==key) return ib->data;
  return NULL;
}

int64_t ih_getint64(struct inthash *ih, int64_t key) {
  struct intbucket *ib = _ih_getval(ih, key);
  if (ib->key==key) return (int64_t)ib->data;
  return IH_INVALID;
}

void ih_setint64(struct inthash *ih, int64_t key, int64_t value) {
  ih_setval(ih,key,(void *)value);
}

void ih_setval(struct inthash *ih, int64_t key, void *data) {
  struct intbucket *ib;
  ib = _ih_getval_deleted(ih, key);
  if (ib->key>IH_SKIP) {
    if (ih->elems>=ih->num_buckets*MAX_LOAD_FACTOR) {
      _ih_add_more_buckets(ih,1);
      ib = _ih_getval(ih, key);
    }
    ih->elems++;
    ib->key = key;
  }
  ib->data = data;
}

void ih_delval(struct inthash *ih, int64_t key) {
  struct intbucket *ib = _ih_getval(ih, key);
  if (key > IH_SKIP || ib->key != key) return;
  ib->key = IH_DELETED;
  ih->elems--;
}

void free_inthash(struct inthash *ih) {
  if (!ih) return;
  if (ih->buckets) {
    memset(ih->buckets, 0, sizeof(struct intbucket)*(ih->num_buckets));
    free(ih->buckets);
  }
  memset(ih, 0, sizeof(struct inthash));
  free(ih);
}

void free_inthash2(struct inthash *ih) {
  int64_t i;
  if (!ih) return;
  for (i=0; i<ih->num_buckets; i++)
    if (ih->buckets[i].key) free_inthash(ih->buckets[i].data);
  free_inthash(ih);
}


void ih_setval2(struct inthash *ih, int64_t key1, int64_t key2, void *data) {
  struct inthash *ih2 = ih_getval(ih, key1);
  if (!ih2) {
    ih2 = new_inthash();
    ih_setval(ih, key1, ih2);
  }
  ih_setval(ih2, key2, data);
}

void *ih_getval2(struct inthash *ih, int64_t key1, int64_t key2) {

  struct inthash *ih2 = ih_getval(ih, key1);
  if (!ih2) return NULL;
  return ih_getval(ih2, key2);
}
