#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ADT.h"

extern int P;

/* Initializes vectorHash structure with k-randomly selected vectors v, k-randomly selected integers r and an integer t in [0,w) */
void createVectorHash(vectorHash *vHash, int k, int dimension)
{
    int i, j;
    
    /* Allocating memory and randomizing vectors v in [0,1] */
    vHash->v = malloc(k*sizeof(int*));
    if (vHash->v == NULL)
    {
        perror("ADT.c: createVectorHash: malloc failed");
        exit(-1);
    }
    for (i = 0; i < k; i++)
    {
        vHash->v[i] = malloc(dimension*sizeof(int));
        if (vHash->v[i] == NULL)
        {
            perror("ADT.c: createVectorHash: malloc failed");
            exit(-1);
        }
        
        for (j = 0; j < dimension; j++)
        {
            vHash->v[i][j] = rand()%2;          /* Vector can only take two values: 0 or 1 */
        }
    }
    
    /* Allocating memory and initializing variables ri */
    vHash->r = malloc(k*sizeof(int));
    if (vHash->r == NULL)
    {
        perror("ADT.c: createVectorHash: malloc failed");
        exit(-1);
    }
    for (i = 0; i < k; i++)
        vHash->r[i] = rand()%10000;
    
    /* t is element of [0,w) ~ w has default value of 4 */
    vHash->t = rand()%W;
}

/* Deallocates memory from vectorHash structure */
void destroyVectorHash(vectorHash *vHash, int k)
{
    int i;
    
    for (i = 0; i < k; i++)
    {
        free(vHash->v[i]);
        vHash->v[i] = NULL;
    }
    
    free(vHash->v);
    vHash->v = NULL;
    free(vHash->r);
    vHash->r = NULL;
}