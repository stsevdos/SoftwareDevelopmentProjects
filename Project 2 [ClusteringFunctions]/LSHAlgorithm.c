/* LSH Algorithm.c */
#include "LSHAlgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "ADT.h"
#include "LSHfunctions.h"

int** LSHAlgorithm(int dimension, int itemCount, char fileType, int k, int L, double **itemsIndex)
{
    item temp;
    int **kDiffNums, i, m, N, j, *hFunctions, *disc, **bucketIndex;
    double ***hx1x2;
    vectorHash* vHash;
    dMatrixItem dmi;
    dMatrixHash* dmh;
       
    N = dimension;

    /**************************************************************************/
    /**** Initialization of structures used for Locality Sensitive Hashing ****/
    if (fileType == 'h')
    {
        /* Create L arrays of k-randomly selected bits */
        kDiffNums = malloc(L*sizeof(int*));
        if (kDiffNums == NULL)
        {
            printf("LSHAlgorithm.c: failed to allocate memory for kDiffNums\n");
            exit(-1);
        }
        /* For each of the L functions g, select the k-bits of hamming item to create each hash function g */
        for(i = 0; i < L; i++)
        {
            kDiffNums[i] = selectKBits(N, k);
        }
    }
    else if (fileType == 'v')
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
            createVectorHash(&(vHash[i]), k, N);
        }
    }
    else /* distance matrix */
    {
    	dmi.distances = itemsIndex;
        /* Create L dMatrixHash structures, each holds data associated with only 1 hash table */
        dmh = malloc(L*sizeof(dMatrixHash));
        if (dmh == NULL)
        {
            printf("LSHAlgorithm.c: failed to allocate memory for dmh\n");
            exit(-1);
        }
        
        /* Create L hx1x2 2d arrays */
        hx1x2 = malloc(L*sizeof(double**));
        if (hx1x2 == NULL)
        {
            printf("LSHAlgorithm: failed to allocate memory for hx1x2\n");
            exit(-1);
        }
        
        /* For L hash functions g */
        for (i = 0; i < L; i++)
        {  
            createDMatrixHash(&(dmh[i]), k, N);          /* Creates structure and selects k random x1, x2 items */
			hx1x2[i] = dbhFucntionH(dmh[i].x1, dmh[i].x2, dmi.distances, N, k);      /* Computes k-h[x1x2] for all items */
			selectt1t2(&dmh[i], hx1x2[i], N, k);          /* Randomly selects k pairs of t1, t2 */
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

    temp.indexNo = -1;
    i = 0;
    for (m = 0; m < itemCount; m++)
    {
        if (fileType == 'h')
        {
            /* Secondary index for item temp structure to be used in hash functions */
            temp.itemCoordinates = itemsIndex[m];
            
            temp.indexNo++;
            /* Hash it through L hashing functions */
            for (i = 0; i < L; i++)
            {
                /* Item #(m) */
                bucketIndex[i][m] = hashFunctionG(temp, kDiffNums[i], k);           /* Get hash table bucket number/index */
            }
        }
        else if (fileType == 'v')
        {
            /* Secondary index for item temp structure to be used in hash functions */
            temp.itemCoordinates = itemsIndex[m];
            
            temp.indexNo++;
            /* Hash it through L hashing functions */
            for (i = 0; i < L; i++)
            {
                hFunctions = euclideanHashFunctionH(vHash[i], temp, N, k);                  /* Pass item through k h functions */
                bucketIndex[i][m] = euclideanHashFunctionF(vHash[i], k, hFunctions);        /* Get hash table bucket number/index */
            }
        }
        else /* distance matrix */
        {
            indexItems(&temp, itemsIndex, NULL, m);                     /* Index each item */
            for (j = 0; j < L; j++)
            {
                disc = discritizeH(dmh[j], hx1x2[j], N, k, i);              /* Discritize each of the k h functions */
                bucketIndex[j][m] = DBHhashFunctionG(disc, k);                      /* Get hash table bucket number/index */
            }
        }
    }
    
    /* Deallocate memory used for any of the 3 files types */
    if (fileType == 'h')      /* Deallocate memory used for Hamming space hashing parameters */
    {
        for (i = 0; i < L; i++)
        {
            free(kDiffNums[i]);
        }
        free(kDiffNums);
    }
    else if (fileType == 'v')     /* Deallocate memory used for Euclidean space hashing parameters */
    {
        free(hFunctions);
        for (i = 0; i < L; i++)
            destroyVectorHash(&vHash[i], k);
        
        free(vHash);
    }
    else            /* Deallocate memory used for DBH parameters */
    {
        destroyDMatrixHash(dmh, k);
        //free(disc);
        
        for (i = 0; i < L; i++)
        {
            for (j = 0; j < N; j++)
                free(hx1x2[i][j]);
            
            free(hx1x2[i]);
        }
        free(hx1x2);
    }
    return bucketIndex;
}
