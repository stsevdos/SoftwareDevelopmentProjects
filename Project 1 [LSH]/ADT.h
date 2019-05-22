/* Abstract Data Types.h */
#ifndef LSH_ADT_h
#define LSH_ADT_h

typedef struct item
{
    double *itemCoordinates;
    char itemName[20];
    int indexNo;
}item;

typedef struct listNode* listNodePtr;
typedef struct listNode
{
    item hashedItem;
    listNodePtr next;
}listNode;

typedef struct vectorHash
{
	int **v;                /* k-randomly selected vectors v for hash function h(p) */
	unsigned int *r;		/* k-randomly selected integers r for hash function f(p) */
	int t;                  /* Randomly selected integer t for hash function h(p) */
}vectorHash;

typedef struct dMatrixItem
{
    double **distances;
    char **dmItemName;
}dMatrixItem;

typedef struct dMatrixHash
{
    int *x1, *x2;       /* k-randomly selected x1, x2 items from distance matrix */
    double *t1, *t2;      /* 2 randomly selected t1, t2 thresholds */
    
}dMatrixHash;

void createList(listNodePtr *start);
void destroyList(listNodePtr *start);
void insertFirst(listNodePtr *start, item input);

listNodePtr** createHashTables(int k, int L);

void createVectorHash(vectorHash *vHash, int k, int dimension);
void destroyVectorHash(vectorHash *vHash, int k);

void createDMatrixHash(dMatrixHash *dMHash, int k, int dimension);
void destroyDMatrixHash(dMatrixHash *dMHash, int k);

#endif