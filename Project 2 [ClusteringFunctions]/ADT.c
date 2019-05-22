/* Abstract Data Types.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ADT.h"

#define W 4

void createList(listNodePtr *start)
{
    *start = NULL;
}

void destroyList(listNodePtr *start)
{
    listNodePtr temp;
    
    while (*start != NULL)
    {
        temp = *start;
        *start = temp->next;
        free(temp);
    }
    *start = NULL;
}

void insertFirst(listNodePtr *start, item input)
{
    listNodePtr newNode = malloc(sizeof(listNode));
    if (newNode == NULL)
    {
        printf("ADT.c: insertFirst: failed to allocate memory for newNode\n");
        exit(-1);
    }
    
    memcpy(&(newNode->hashedItem), &input, sizeof(item));
    newNode->next = *start;
    *start = newNode;
}

/* Function that creates a table of L hash functions g */
listNodePtr** createHashTables(int k, int L)
{
    int i, j, hashTableSize;
    listNodePtr** ht;
    
    hashTableSize = (int)(0.5 + pow(2, k)); /* in case we get e.g. 8, it's 7.99999 and gets trancated to 7 */
    
    /* Allocates memory for L hash tables and initializes each of L indexes to list pointer */
    ht = malloc(L*sizeof(listNodePtr*));
    if (ht == NULL)
    {
        printf("ADTc.: createHashTables: failed to allocate memory for ht\n");
        exit(-1);
    }
    
    for (i = 0; i < L; i++)
    {
        ht[i] = malloc(hashTableSize * sizeof(listNodePtr));
        if (ht[i] == NULL)
        {
            printf("ADT.c: createHashTables: failed to allocate memory for ht[%d]\n)", i);
            exit(-1);
        }
        
        for (j = 0; j < hashTableSize; j++)
        {
            createList(&(ht[i][j]));
        }
    }
    return ht;
}

/* Initializes vectorHash structure with k-randomly selected vectors v, k-randomly selected integers r and an integer t in [0,w) */
void createVectorHash(vectorHash *vHash, int k, int dimension)
{
    int i, j;
    
    /* Allocating memory and randomizing vectors v in [0,1] */
    vHash->v = malloc(k*sizeof(int*));
    if (vHash->v == NULL)
    {
        printf("ADT.c: createVectorHash: failed to allocate memory for vHash->v\n");
        exit(-1);
    }
    for (i = 0; i < k; i++)
    {
        vHash->v[i] = malloc(dimension*sizeof(int));
        if (vHash->v[i] == NULL)
        {
            printf("ADT.c: createVectorHash: failed to allocate memory for vHash->[%d]\n", i);
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
        printf("ADT.c: createVectorHash: failed to allocate memory for vHash->r\n");
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

/* Initializes struct with hashing information except for thresholds t1, t2 */
void createDMatrixHash(dMatrixHash *dMHash, int k, int dimension)
{
    int i, j, sameExists;
    
    /* Allocates memory for arrays for x1 and x2 random items p (their array indexes) */
    dMHash->x1 = malloc(k*sizeof(int));
    if (dMHash->x1 == NULL)
    {
        printf("ADT.c: createDMatrixHash: failed to allocate memory for dMHash->x1\n");
        exit(-1);
    }
    dMHash->x2 = malloc(k*sizeof(int));
    if (dMHash->x2 == NULL)
    {
        printf("ADT.c: createDMatrixHash: failed to allocate memory for dMHash->x2\n");
        exit(-1);
    }
    
    /*  Randomly selects k-pairs of x1, x2 items */
    for (i = 0; i < k;)
    {
        dMHash->x1[i] = rand()%dimension;		/* Holds index of item */
        do
        {
            dMHash->x2[i] = rand()%dimension;	/* Holds index of item */
        }while(dMHash->x1[i] == dMHash->x2[i]);
        
        /* Checks for equality between pairs of x1, x2 */
        sameExists = 0;
        for (j = 0; j < i; j++)
        {
            if (((dMHash->x1[j] == dMHash->x1[i]) && (dMHash->x2[j] == dMHash->x2[i])) || ((dMHash->x1[j] != dMHash->x1[i]) && (dMHash->x2[j] != dMHash->x2[i])))
                sameExists = 1;
        }
        if (sameExists == 0)
            i++;
    }
}

/* Deallocates memoery allocated for dMatrixHash structure */
void destroyDMatrixHash(dMatrixHash *dMHash, int k)
{
    free(dMHash->x1);
    dMHash->x1 = NULL;
    free(dMHash->x2);
    dMHash->x2 = NULL;
    free(dMHash->t1);
    dMHash->t1 = NULL;
    free(dMHash->t2);
    dMHash->t2 = NULL;
}