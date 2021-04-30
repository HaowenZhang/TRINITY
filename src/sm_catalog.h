#ifndef SM_CATALOG_H
#define SM_CATALOG_H

#define ND_START -8
#define ND_STOP 1
#define ND_BPDEX 10
#define ND_INV_BPDEX (1.0/(double)ND_BPDEX)
#define ND_BINS ((ND_STOP-ND_START)*ND_BPDEX + 2)

#define Z_START 0
#define Z_STOP 12
#define Z_BPZ 10
#define Z_INV_BPZ (1.0/(double)Z_BPZ)
#define Z_BINS ((Z_STOP-Z_START)*Z_BPZ + 2)

#endif /* SM_CATALOG_H */
