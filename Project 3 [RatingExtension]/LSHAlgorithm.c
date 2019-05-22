/* LSH Algorithm.c */
#include "LSHAlgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "ADT.h"
#include "LSHFunctions.h"

extern int P;

short** LSHAlgorithm(int dimension, int itemCount, char fileType, int k, int L, short **ratings, float *r)
{
    int i, m, j, temp, *hFunctions, **bucketIndex;
    short ***csVectors, neighbourCounter, **pNearest;
    double radius, initRadius, vicinity;
    vectorHash* vHash;
    
    /**************************************************************************/
    /**** Initialization of structures used for Locality Sensitive Hashing ****/
    if (fileType == 'v')
    {
        /* Allocate memory for L vector hash structures, containing randomly selected v vectors in [0,1], t randomly selected in [0, w) and randomly selected ri integers for hash function g */
        vHash = malloc(L*sizeof(vectorHash));
        if (vHash == NULL)
        {
            printf("LSHAlgorithm.c: failed to allocate memory for vHash\n");
            exit(-1);
        }
        /* Creating L vector hash structures */
        for (i = 0; i < L; i++)
        {
            createVectorHash(&(vHash[i]), k, dimension);
        }
        initRadius = 3.74;
        radius = initRadius;
    }
    else /* Cosine similarity */
    {
        csVectors = malloc(sizeof(short**)*L);
        /* check */
        for (i = 0; i < L; i++)
        {
            csVectors[i] = cosineRandomVectors(k, dimension);
        }
    }
    
    bucketIndex = malloc(sizeof(int*)*L);
    if (bucketIndex == NULL)
    {
        printf("LSHAlgorithm: failed to allocate memory for bucketIndex\n");
        exit(-1);
    }
    for (i = 0; i < L; i++)
    {
        bucketIndex[i] = malloc(sizeof(int)*itemCount);
        if (bucketIndex[i] == NULL)
        {
            printf("LSHAlgorithm: failed to allocate memory for bucketIndex[%d]\n", i);
            exit(-1);
        }
    }
    /**** End Initialization of structures used for Locality Sensitive Hashing ****/
    /*****************************************************************************/
    
    /********** Creating Hash Tables Index **********/
    for (m = 0; m < itemCount; m++)
    {
        if (fileType == 'v')
        {
            /* Hash it through L hashing functions */
            for (i = 0; i < L; i++)
            {
                hFunctions = euclideanHashFunctionH(vHash[i], ratings[i], dimension, k);                  /* Pass item through k h functions */
                bucketIndex[i][m] = euclideanHashFunctionF(vHash[i], k, hFunctions);        /* Get hash table bucket number/index */
            }
        }
        else /* Cosine Similarity */
        {
            for (i = 0; i < L; i++)
            {
                csVectors[i] = cosineRandomVectors(k, dimension);
                bucketIndex[i][m] = cosineHashFunction(k, dimension, csVectors[i], ratings[m], r);
            }
        }
    }
    /********** Finding P nearest neighbours **********/
    pNearest = malloc(sizeof(short*)*itemCount);    /* check */
    
    for (i = 0; i < itemCount; i++)
    {
        pNearest[i] = malloc(sizeof(short)*P);  /* check */
        for (j = 0; j < P; j++)
            pNearest[i][j] = -1.0;
        
        radius = initRadius;
        for (neighbourCounter = 0; neighbourCounter < P && radius <= 1.0; radius *= 1.5)
        {
            for (j = 0; j < L; j++)
            {
                for (m = 0; m < itemCount; m++)
                {
                    if (m == i)
                        continue;
                    /* If item m is in same bucket as item i */
                    if (bucketIndex[j][i] == bucketIndex[j][m])
                    {
                        if (fileType == 'v')
                            vicinity = euclideanDistance(i, m, ratings, dimension);
                        else
                            vicinity = cosineSimilarity(i, m, ratings, dimension);
                        
                        if (vicinity <= radius)
                        {
                            for (temp = 0; temp < neighbourCounter; temp++)
                            {
                                if (pNearest[i][temp] == m)
                                    break;
                            }
                            if (temp < neighbourCounter)
                                continue;
                            
                            pNearest[i][neighbourCounter] = m;
                            neighbourCounter++;
                        }
                    }
                }
            }
        }
    }
   
    /* Deallocate memory used for any of the 3 files types */
    if (fileType == 'v')     /* Deallocate memory used for Euclidean space hashing parameters */
    {
        free(hFunctions);
        for (i = 0; i < L; i++)
            destroyVectorHash(&vHash[i], k);
        
        free(vHash);
    }
    else            /* Deallocate memory used for cosine similarity parameters */
    {
        for (i = 0; i < L; i++)
        {
            for (j = 0; j < k; j++)
                free(csVectors[i][j]);
            free(csVectors[i]);
        }
    }
    
    for (i = 0; i < L; i++)
        free(bucketIndex[i]);
    free(bucketIndex);
    
    return pNearest;
}
