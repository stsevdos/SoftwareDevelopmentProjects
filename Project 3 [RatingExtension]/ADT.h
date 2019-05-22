#ifndef __Recommendation__ADT__
#define __Recommendation__ADT__

#include <stdio.h>

#define W 4

typedef struct vectorHash
{
    int **v;                /* k-randomly selected vectors v for hash function h(p) */
    unsigned int *r;		/* k-randomly selected integers r for hash function f(p) */
    int t;                  /* Randomly selected integer t for hash function h(p) */
}vectorHash;

void createVectorHash(vectorHash *vHash, int k, int dimension);
void destroyVectorHash(vectorHash *vHash, int k);

#endif /* defined(__Recommendation__ADT__) */
