/* Functions for Locality Sensitive Hashing */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ADT.h"
#include "LSHfunctions.h"

#define W 4

/* Randomly selects a sequence of k-bits to define hamming hashing functions ~ N = dimension, k = random bits */
int* selectKBits(int N, int k)
{
    int i, j, sameExists;
    int *indexes = malloc(k*sizeof(int));
    if (indexes == NULL)
    {
        printf("LSHfunctions.c: selectKBits: failed to allocate memory for indexes\n");
        exit(-1);
    }
    
    indexes[0] = rand() % N;            /* Initializes first bit outside of loop so it can compare two of the for equality inside loop */
    for (i = 1; i < k;)
    {
        indexes[i] = rand() % N;
        
        sameExists = 0;
        for (j = 0; j < i; j++)         /* Checking if array indexes selected up to i have been selected before and if so randomly selects again */
        {
            if (indexes[i] == indexes[j])
                sameExists = 1;
        }
        if (sameExists == 0)           /* if same array index doesn't exist in already randomly selected k-bts, continue with next random bit */
            i++;
    }
    return indexes;                     /* returns the array of k-randomly selected bits (in array indexes, not actual hamming bits) to hash point in bucket */
}

/* Returns index (bucket) for input item to be hashed in */
int hashFunctionG(item input, int* kBits, int k)
{
    int i, itemIndex = 0;
    
    /* Computes selected bits from binary to decimal to determine bucket (array index) */
    for (i = 0; i < k; i++)
        itemIndex = (itemIndex * 2) + (int)input.itemCoordinates[kBits[i]];
    
    return itemIndex;
}

/* Reads 1 input item from file (hamming file), N = dimension, returns 0 end of file reached, 1 otherwise */
int readHamming(FILE* inputFile, int N, item* currentItem)
{
    char i, itemCoord[70];
    
    /* If reached end of file, return 0 */
    if (feof(inputFile))
        return 0;
    
    /* Reading item name */
    fscanf(inputFile, "%s", currentItem->itemName);
    /* Reading item coordinates */
    fscanf(inputFile, "%s", itemCoord);
    
    /* Hamming coordinates are read as a string of doubles */
    currentItem->itemCoordinates = malloc(N*sizeof(double));
    if (currentItem->itemCoordinates == NULL)
    {
        printf("LSHfunctions.c: readHamming: failed to allocate memory for currentItem->itemCoordinates\n");
        exit(-1);
    }
    
    /* Copies bits read to item structure */
    for(i = 0; i < N; i++)
        currentItem->itemCoordinates[i] = (double) (itemCoord[i] - '0');
    
    return 1;
}

/* Reads 1 input item from file (euclidean/vector file), N = dimension, returns 0 end of file reached, 1 otherwise */
int readVector(FILE* inputFile, int N, item* currentItem)
{
    int i;
    
    if (feof(inputFile))	/* If reached end of file, return 0 */
        return 0;
    
    fscanf(inputFile, "%s", currentItem->itemName);         /* Reading item name */
    currentItem->itemCoordinates = malloc(N*sizeof(double));
    if (currentItem->itemCoordinates == NULL)
        printf("LSHfunctions.c: readVector: failed to allocate memory for currentItem->itemCoordinates\n");
    
    /* Reading item coordinates */
    for(i = 0; i < N; i++)
    {
        fscanf(inputFile, "%lf", &(currentItem->itemCoordinates[i] ));
    }
    
    return 1;
}


/* Reads first item from file and returns its dimension */
int getHammingDimension(char* fileName)
{
    FILE *inputFile;
    int dimension;
    char* stringBuffer;
    
    stringBuffer = malloc(sizeof(char)*100);        /* max number of bits should be 64 */
    if (stringBuffer == NULL)
    {
        printf("LSHfunctions.c: getHammingDimension: failed to allocate memory for stringBuffer\n");
        exit(-1);
    }
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        printf("LSHfunctions.c: getHammingDimension: failed to open file (inputFile) %s\n", fileName);
        exit(-1);
    }
    
    /* Read first line from file */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    /* Read item name */
    fscanf(inputFile, "%s", stringBuffer);
    /* Read item coordinates */
    fscanf(inputFile, "%s", stringBuffer);
    
    dimension = strlen(stringBuffer);       /* max should be 64 */
    
    fclose(inputFile);
    free(stringBuffer);
    stringBuffer = NULL;
    
    return dimension;
}

/* Finds the distance between two given items in hamming space */
int hammingDistance(item input, item query, int dimension)
{
    int distance = 0, i;
    
    for (i = 0; i < dimension; i++)
    {
        /* For each pair of bits that are different, increment the distance value */
        if (input.itemCoordinates[i] != query.itemCoordinates[i])
            distance++;
    }
    
    return distance;
}

/* Returns the dimension of a vector */
int getVectorDimension(char* fileName, char* metric)
{
    FILE *inputFile;
    
    int dimension = 1;      /* Initialize to 1 because we miss the first blank space before 1st coordinate */
    char charBuffer;
    char* stringBuffer;
    
    stringBuffer = malloc(sizeof(char)*100);
    if (stringBuffer == NULL)
    {
        printf("LSHfunctions.c: getVectorDimension: failed to allocate memory for stringBuffer\n");
        exit(-1);
    }
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        printf("LSHfunctions.c: getVectorDimension: failed to open file (inputFile) %s\n", fileName);
        exit(-1);
    }
    
    /* Read first line from file */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    /* Read second line from file */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    *metric = stringBuffer[0];
    /* Read item name */
    fscanf(inputFile, "%s", stringBuffer);
    
    /* Each blank space read increments dimension */
    fscanf(inputFile, "%c", &charBuffer);
    while (charBuffer != '\n')
    {
        fscanf(inputFile, "%c", &charBuffer);
        if (charBuffer == ' ' || charBuffer == '\t')
            dimension++;
    }
    
    fclose(inputFile);
    free(stringBuffer);
    stringBuffer = NULL;
    
    return dimension;
}

/* Returns array of all h(p) functions for a single point in space */
int* euclideanHashFunctionH(vectorHash vh, item vectorP, int dimension, int k)
{
    int i, j;
    
    int *hFunctions = malloc(sizeof(int)*k);
    if (hFunctions == NULL)
    {
        printf("LSHfunctions.c: euclideanHashFunctionH: failed to allocate memory for hFunctions\n");
        exit(-1);
    }
    
    /* Calculating k-hi functions for one g hash function */
    for (i = 0; i < k; i++)
    {
        hFunctions[i] = 0;
        for (j = 0; j < dimension; j++)
        {
            if (vh.v[i][j] == 1)
                hFunctions[i] += vectorP.itemCoordinates[j];
        }
        
        hFunctions[i] = (hFunctions[i] + vh.t)/W;
    }
    
    return hFunctions;
}

/* Returns index of hash table for a single point in space */
int euclideanHashFunctionF(vectorHash vh, int k, int *hFunctions)
{
    int i, f = 0, tableSize;
    
    long M;
    
    M = pow(2, 31) - 5;
    
    tableSize = 0.5 + pow(2, k);
    
    for (i = 0; i < k; i++)
    {
        f += (hFunctions[i]*vh.r[i]) % M;
    }
    
    f = f % tableSize;
    
    return f;
}

/* Calculates and returns euclidean distance between two vectors in the euclidean space */
double euclideanDistance(item input, item query, int dimension)
{
    double distance = 0;
    int i;
    
    for (i = 0; i < dimension; i++)
    {
        distance += (input.itemCoordinates[i] - query.itemCoordinates[i])*(input.itemCoordinates[i] - query.itemCoordinates[i]);
    }
    
    distance = sqrt(distance);
    
    return distance;
}

/* Returns dimension of a given distance matrix */
int getDistanceMatrixDimension(char* fileName)
{
    FILE *inputFile;
    
    int dimension = 0;
    char charBuffer;
    char* stringBuffer;
    
    stringBuffer = malloc(sizeof(char)*100);
    if (stringBuffer == NULL)
    {
        printf("LSHfunctions.c: getDistanceMatrixDimension: failed to allocate memory for stringBuffer\n");
        exit(-1);
    }
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        printf("LSHfunctions.c: getDistanceMatrixDimension: failed to open file (inputFile) %s\n", fileName);
        exit(-1);
    }
    
    /* Read first line from file */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    
    /* Read second line from file */
    while (charBuffer != '\n')
        fscanf(inputFile, "%c", &charBuffer);
    
    charBuffer = 'a';       /* charBuffer previously held \n, so it would never enter the while loop below */
    while (charBuffer != '\n')
    {
        /* Each blank space or tab character met increments dimension value */
        fscanf(inputFile, "%c", &charBuffer);
        if (charBuffer == ' ' || charBuffer == '\t')
            dimension++;
    }
    
    fclose(inputFile);
    free(stringBuffer);
    stringBuffer = NULL;
    
    return dimension;
}

/* Reads item names and distances from file and inserts them in corresponding matrices */
void readDistances(FILE* inputFile, int dimension, dMatrixItem *dmi)
{
    int i, j;
    char *stringBuffer;
    
    stringBuffer = malloc(sizeof(char)*30);
    if (stringBuffer == NULL)
    {
        printf("LSHfunctions.c: readDistances: failed to allocate memory for stringBuffer\n");
        exit(-1);
    }
    
    dmi->dmItemName = malloc(sizeof(char*)*dimension);
    if (dmi->dmItemName == NULL)
    {
        printf("LSHfunctions.c: readDistances: failed to allocate memory for dmi->dmItemName\n");
        exit(-1);
    }
    
    /* Reading first string from second line */
    fscanf(inputFile, "%s", stringBuffer);
    
    /* Reading items' names and inserting them into an array of strings */
    for (i = 0; i < dimension; i++)
    {
        dmi->dmItemName[i] = malloc(sizeof(char)*30);
        if (dmi->dmItemName[i] == NULL)
        {
            printf("LSHfunctions.c: readDistances: failed to allocate memory for dmi->dmItemName[%d]\n", i);
            exit(-1);
        }
        fscanf(inputFile, "%s", dmi->dmItemName[i]);
    }
    
    /* Reading distances from file and inserting them into 2-d array */
    dmi->distances = malloc(sizeof(double*)*(dimension+1));
    if (dmi->distances == NULL)
    {
        printf("LSHfunctions.c: readDistances: failed to allocate memory for dmi->distances\n");
        exit(-1);
    }
    for (i = 0; i < dimension; i++)
    {
        dmi->distances[i] = malloc(sizeof(double)*dimension);
        if (dmi->distances[i] == NULL)
        {
            printf("LSHfunctions.c: readDistances: failed to allocate memory for dmi->distances[%d]\n", i);
            exit(-1);
        }
        for (j = 0; j < dimension; j++)
        {
            fscanf(inputFile, "%lf", &(dmi->distances[i][j]));
        }
    }
    
    dmi->distances[dimension] = malloc(sizeof(double)*dimension);
    if (dmi->distances[dimension] == NULL)
    {
        printf("LSHfunctions.c: readDistances: failed to allocate memory for dmi->distances[%d]\n", dimension);
        exit(-1);
    }
    
    /* in case of getting an upper triangular matrix, copy distances to 0s */
    for (i = 0; i < dimension; i++)
        for (j = 0; j < i; j++)
            dmi->distances[i][j] = dmi->distances[j][i];
    
    free(stringBuffer);
    stringBuffer = NULL;
}

/* Function that indexes the items for hashing */
void indexItems(item *temp, double **itemDistances, char **itemNames, int indexPosition)
{
    temp->itemCoordinates = itemDistances[indexPosition];
    strcpy(temp->itemName, itemNames[indexPosition]);
    temp->indexNo = indexPosition;
}

/* Calculates all functions hi for all points and returns a 2-d array */
double** dbhFucntionH(int *x1, int *x2, double **distancesMatrix, int dimension, int k)
{
    int i, j;
    double **hx1x2;
    
    hx1x2 = malloc(sizeof(double*)*dimension);
    if (hx1x2 == NULL)
    {
        printf("LSHfunctions.c: dbhFucntionH: failed to allocate memory for hx1x2\n");
        exit(-1);
    }
    
    for (i = 0; i < dimension; i++)
    {
        hx1x2[i] = malloc(sizeof(double)*k);
        if (hx1x2[i] == NULL)
        {
            printf("LSHfunctions.c: dbhFunctionH: failed to allocate memory for hx1x2[%d]\n", i);
            exit(-1);
        }
        
        for (j = 0; j < k; j++)
        {
            
            hx1x2[i][j] = (distancesMatrix[i][x1[j]]*distancesMatrix[i][x1[j]] + distancesMatrix[i][x2[j]]*distancesMatrix[i][x2[j]] - distancesMatrix[x1[j]][x2[j]]*distancesMatrix[x1[j]][x2[j]])/(2*distancesMatrix[x1[j]][x2[j]]);
        }
    }
    
    return hx1x2;
}

/* Selection sort implementation */
void SelectionSort(double* array, int arrayLength)
{
    int i, j;
    double temp;
    
    for (i = 0; i < arrayLength - 1; i++)
    {
        for (j = i+1; j < arrayLength; j++)
        {
            if (array[i] > array[j])
            {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }
    }
}

/* Sorts hx1x2 matrix and selects t1, t2 thresholds with t1 = median(hx1x2), t2 = max(hx1x2) */
void selectt1t2(dMatrixHash *dmh, double **hx1x2, int dimension, int k)
{
    int i, j;
    int median;
    double *hp;
    
    hp = malloc(sizeof(double)*dimension);
    if (hp == NULL)
    {
        printf("LSHfunctions.c: selectt1t2: failed to allocate memory for hp\n");
        exit(-1);
    }
    
    /* median of hx1x2 is its dimension/2 */
    median = dimension/2;
    
    dmh->t1 = malloc(sizeof(double)*k);
    if (dmh->t1 == NULL)
    {
        printf("LSHfunctions.c: selectt1t2: failed to allocate memory for dmh->t1\n");
        exit(-1);
    }
    dmh->t2 = malloc(sizeof(double)*k);
    if (dmh->t2 == NULL)
    {
        printf("LSHfunctions.c: selectt1t2: failed to allocate memory for dmh->t2\n");
        exit(-1);
    }
    
    /* for each of the k columns, sort it and then select t1 = median(hx1x2) and t2 = max(hx1x2) */
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            hp[j] = hx1x2[j][i];
        }
        
        SelectionSort(hp, dimension);
        
        dmh->t1[i] = hp[median];
        dmh->t2[i] = hp[dimension-1] +1 ;
    }
    
    free(hp);
    hp = NULL;
}

/* Discritizes functions' (hx1x2) values based on thresholds (t1, t2) and returns array of 0s and 1s */
int* discritizeH(dMatrixHash dmh, double **hx1x2, int dimension, int k, int index)
{
    int i;
    int *discritize;
    
    discritize = malloc(sizeof(int)*k);
    if (discritize == NULL)
    {
        printf("LSHfunctions.c: discritizeH: failed to allocate memory for discritize\n");
        exit(-1);
    }
    
    /* if hx1x2 is element of [t1, t2] then h[(x1x2),(t1,t2)] is 1, otherwise 0 */
    for (i = 0; i < k; i++)
    {
        if (hx1x2[index][i] >= dmh.t1[i] && hx1x2[index][i] <= dmh.t2[i])
            discritize[i] = 1;
        else
            discritize[i] = 0;
    }
    
    return discritize;
}

/* Calculates and returns index position for item to be hashed */
int DBHhashFunctionG(int* kBits, int k)
{
    int i, itemIndex = 0;
    
    for (i = 0; i < k; i++)
        itemIndex = (itemIndex * 2) + kBits[i];
    
    return itemIndex;
}

/* Processes query item hashing and returns bucket (array index) it was hashed */
int queryHashingDBH(double **distancesMatrix, int dimension, dMatrixHash dmh, int k, int hashTableNo)
{
    int j, hashIndex;
    int *discritize;
    double *hx1x2Query;
    
    discritize = malloc(sizeof(int)*k);
    if (discritize == NULL)
    {
        printf("LSHfunctions.c: queryHashingDBH: failed to allocate memory for discritize\n");
        exit(-1);
    }
    
    hx1x2Query = malloc(sizeof(double)*k);
    if (hx1x2Query == NULL)
    {
        printf("LSHfunctions.c: queryHashingDBH: failed to allocate memory for hx1x2Query\n");
        exit(-1);
    }
    
    /* Calculate hx1x2 fucntions for query item */
    for (j = 0; j < k; j++)
    {
        hx1x2Query[j] = (distancesMatrix[dimension][dmh.x1[j]]*distancesMatrix[dimension][dmh.x1[j]] + distancesMatrix[dimension][dmh.x2[j]]*distancesMatrix[dimension][dmh.x2[j]] - distancesMatrix[dmh.x1[j]][dmh.x2[j]]*distancesMatrix[dmh.x1[j]][dmh.x2[j]])/(2*distancesMatrix[dmh.x1[j]][dmh.x2[j]]);
    }
    
    /* Discritize hx1x2 functions to 0 and 1 */
    for (j = 0; j < k; j++)
    {
        if (hx1x2Query[j] >= dmh.t1[j] && hx1x2Query[j] <= dmh.t2[j])
            discritize[j] = 1;
        else
            discritize[j] = 0;
    }
    
    /* Get hash table index for query item */
    hashIndex = DBHhashFunctionG(discritize, k);
    
    free(discritize);
    discritize = NULL;
    free(hx1x2Query);
    hx1x2Query = NULL;
    
    return hashIndex;
}