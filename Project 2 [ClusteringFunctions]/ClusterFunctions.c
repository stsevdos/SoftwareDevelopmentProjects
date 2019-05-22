#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ClusterFunctions.h"
#include "ADT.h"
#include "LSHfunctions.h"

void readConfigFile(char *configFileName, int *k, int *hashK, int *L, int *maxNeighbor, int *numLocal)
{
    int tempRead;
    char *stringBuffer;
    FILE *configFile;
    
    stringBuffer = malloc(sizeof(char)*50);
    if (stringBuffer == NULL)
    {
        printf("ClusterFunctions.c: readConfigFile: Failed to allocate memory for stringBuffer\n");
        exit(-1);
    }
    
    configFile = fopen(configFileName, "r");
    if (configFile == NULL)
    {
        printf("ClusterFunctions.c: readConfigFile: Could not open (configuration) file %s\n", configFileName);
        exit(-1);
    }
    while (!feof(configFile))
    {
        fscanf(configFile, "%s", stringBuffer);
        if (strcmp(stringBuffer, "number_of_clusters:") == 0)
        {
            fscanf(configFile, "%d", &tempRead);
            if (tempRead != 0)
                *k = tempRead;
        }
        else if (strcmp(stringBuffer, "number_of_hash_functions:") == 0)
        {
            fscanf(configFile, "%d", &tempRead);
            if (tempRead != 0)
                *hashK = tempRead;
        }
        else if (strcmp(stringBuffer, "number_of_hash_tables:") == 0)
        {
            fscanf(configFile, "%d", &tempRead);
            if (tempRead != 0)
                *L = tempRead;
        }
        else if (strcmp(stringBuffer, "clarans_set_fraction:") == 0)
        {
            fscanf(configFile, "%d", &tempRead);
            if (tempRead != 0)
                *maxNeighbor = tempRead;
        }
        else if (strcmp(stringBuffer, "clarans_iterations:") == 0)
        {
            fscanf(configFile, "%d", &tempRead);
            if (tempRead != 0)
                *numLocal = tempRead;
        }
        else
        {
            printf("ClusterFunctions.c: readConfigFile: Unsupported data\n");
            exit(-1);
        }
    }
    free(stringBuffer);
    fclose(configFile);
}

int getItemCount(char *fileName, char* fileType)
{
    FILE *inputFile;
    
    int itemCount;
    char charBuffer;
    char* stringBuffer;
    
    stringBuffer = malloc(sizeof(char)*100);
    if (stringBuffer == NULL)
    {
        printf("ClusterFunctions.c: getItemCount: failed to allocate memory for stringBuffer\n");
        exit(-1);
    }
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        printf("ClusterFunctions.c: getItemCount: failed to open file (inputFile) %s\n", fileName);
        exit(-1);
    }
    
    /* Read first line from file */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    
    itemCount = 0;
    
    /* If items in file are in vector space, read second line of file as well */
    if (strcmp(stringBuffer, "vector") == 0 || strcmp(stringBuffer, "euclidean") == 0)
    {
        fscanf(inputFile, "%s", stringBuffer);
        fscanf(inputFile, "%s", stringBuffer);
        *fileType = 'v';
    }
    else if (strcmp(stringBuffer, "hamming") == 0)
        *fileType = 'h';
    else
    {
        fscanf(inputFile, "%c", &charBuffer);
        while (charBuffer != '\n')
            fscanf(inputFile, "%c", &charBuffer);
        *fileType = 'm';
    }
    
    while (!feof(inputFile))
    {
        /* Each new line character met increments itemCount value */
        fscanf(inputFile, "%c", &charBuffer);
        if (charBuffer == '\n')
            itemCount++;
    }
    
    fclose(inputFile);
    free(stringBuffer);
    
    return itemCount;
}

double** readItems(char *fileName, char **itemNamesIndex, int dimension, int itemCount)
{
    int i, j;
    double **itemsIndex;
    char *stringBuffer;
    dMatrixItem dmi;
    FILE *inputFile;
    
    stringBuffer = malloc(sizeof(char)*100);
    if (stringBuffer == NULL)
    {
        printf("ClusterFunctions.c: readItems: Failed to allocate memory for stringBuffer\n");
        exit(-1);
    }

    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        printf("ClusterFunctions.c: readItems: Failed to open file (inputFile) %s\n", fileName);
        exit(-1);
    }
    itemsIndex = malloc(sizeof(double*)*itemCount);
    if (itemsIndex == NULL)
    {
        printf("ClusterFunctions.c: readItems: Failed to allocate memory for itemsIndex\n");
        exit(-1);
    }
    for (i = 0; i < itemCount; i++)
    {
        itemsIndex[i] = malloc(sizeof(double)*dimension);
        if (itemsIndex[i] == NULL)
        {
            printf("ClusterFunctions.c: readItems: Failed to allocate memory for itemsIndex[%d]\n", i);
            exit(-1);
        }
        
    }
    /* Read first line from file */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    
    if (strcmp(stringBuffer, "vector") == 0)
    {
        /* If items in file are in vector space, read second line of file as well */
        fscanf(inputFile, "%s", stringBuffer);
        fscanf(inputFile, "%s", stringBuffer);
        
        for (i = 0; i < itemCount; i++)
        {
            /* Read item's name, then its coordinates */
            fscanf(inputFile, "%s", stringBuffer);
            strcpy(itemNamesIndex[i], stringBuffer);
            for (j = 0; j < dimension; j++)
                fscanf(inputFile, "%lf", &itemsIndex[i][j]);
        }
    }
    else if (strcmp(stringBuffer, "hamming") == 0)
    {
        for (i = 0; i < itemCount; i++)
        {
            /* Read item's name, then its coordinates */
            fscanf(inputFile, "%s", stringBuffer);
            strcpy(itemNamesIndex[i], stringBuffer);
            /* Read item's "coordinates" */
            fscanf(inputFile, "%s", stringBuffer);
            /* Copy "coordinates" into itemsIndex */
            for (j = 0; j < dimension; j++)
                itemsIndex[i][j] = (double) (stringBuffer[j] - '0');
        }
    }
    else
    {
        /* readDistances puts file's item coordinates and names into dmi struct */
        readDistances(inputFile, dimension, &dmi);
        itemsIndex = dmi.distances;
        for (i = 0; i < itemCount; i++)
            itemNamesIndex[i] = dmi.dmItemName[i];
    }
    
    free(stringBuffer);
    fclose(inputFile);
    
    return itemsIndex;
}

int getDimension(char *fileName, char fileType)
{
    int dimension;
    
    if (fileType == 'h')
        dimension = getHammingDimension(fileName);
    else if (fileType == 'v')
        dimension = getVectorDimension(fileName, &fileType);
    else
        dimension = getDistanceMatrixDimension(fileName);
    
    return dimension;
}

double** calcDistanceMatrix(double **itemsIndex, char fileType, int dimension, int itemCount)
{
    int i, j;
    double **distanceMatrix;
    char *stringBuffer;
    item temp1, temp2;
    
    stringBuffer = malloc(sizeof(char)*100);
    if (stringBuffer == NULL)
    {
        printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for stringBufffer\n");
        exit(-1);
    }
    distanceMatrix = malloc(sizeof(double*)*itemCount);
    if (distanceMatrix == NULL)
    {
        printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for distanceMatrix\n");
        exit(-1);
    }
    
    if (fileType == 'h')
    {
        /* If items are in Hamming Space */
        for (i = 0; i < itemCount; i++)
        {
            distanceMatrix[i] = malloc(sizeof(double)*itemCount);
            if (distanceMatrix[i] == NULL)
            {
                printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for distanceMatrix[%d]\n", i);
                exit(-1);
            }
            temp1.itemCoordinates = itemsIndex[i];
            /* Use hammingDistance to calculate distance for distance matrix */
            for (j = 0; j < itemCount; j++)
            {
                temp2.itemCoordinates = itemsIndex[j];
                distanceMatrix[i][j] = hammingDistance(temp1, temp2, dimension);
            }
        }
    }
    else if (fileType == 'v')
    {
        /* If items are in Vector Space */
        for (i = 0; i < itemCount; i++)
        {
            distanceMatrix[i] = malloc(sizeof(double)*itemCount);
            if (distanceMatrix[i] == NULL)
            {
                printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for distanceMatrix[%d]\n", i);
                exit(-1);
            }
            temp1.itemCoordinates = itemsIndex[i];
            /* Use euclideanDistnace to calculate distance for distance matrix */
            for (j = 0; j < itemCount; j++)
            {
                temp2.itemCoordinates = itemsIndex[j];
                distanceMatrix[i][j] = (double)euclideanDistance(temp1, temp2, dimension);
            }
        }
    }
    else    /* If items are part of a distance matrix, return itemsIndex pointer itself as distanceMatrix */
        return itemsIndex;
    
    free(stringBuffer);
    
    return distanceMatrix;
}

int sameMedoids(int *oldMedoids, int *newMedoids, int k)
{
    int i, j, notExists;
    
    notExists = 1;
    /* For each old medoid */
    for (i = 0; i < k; i++)
    {
        notExists = 1;
        /* Check if it has been selected before in any position of newMedoids */
        for (j = 0; j < k; j++)
        {
            if (oldMedoids[i] == newMedoids[j])
            {
                /* If found, continue checking for other old medoids */
                notExists = 0;
                break;
            }
        }
        /* If not found, arrays don't have the same elements, return 0 */
        if (notExists == 1)
            return 0;
    }
    return 1;
}

int* initializationConcentrate(double **distanceMatrix, int itemCount, int k)
{
    int i, j;
    int *initMedoids, *vIndexes;
    double *v, *distSums;

    v = malloc(sizeof(double)*itemCount);
    if (v == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for v array\n");
        exit(-1);
    }
    vIndexes = malloc(sizeof(int)*itemCount);
    if (vIndexes == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for vIndexes array\n");
        exit(-1);
    }
    distSums = malloc(sizeof(double)*itemCount);
    if (distSums == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for distSums array\n");
        exit(-1);
    }
    initMedoids = malloc(sizeof(int)*k);
    if (initMedoids == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for initMedoids\n");
        exit(-1);
    }
    for (i = 0; i < itemCount; i++)
    {
        /* Calculating denominator for vi values (for each item) */
        distSums[i] = 0;
        for (j = 0; j < itemCount; j++)
            distSums[i] += distanceMatrix[i][j];
    }
    for (i = 0; i < itemCount; i++)
    {
        v[i] = 0;
        vIndexes[i] = i;    /* Keep item index for sorting */
        /* Calculating vi values for each item */
        for (j = 0; j < itemCount; j++)
            v[i] += (double)(distanceMatrix[i][j]/distSums[j]);
    }
    /* Sort both arrays to insure valid vIndexes values */
    SelectionSortTwoArrays(v, vIndexes, itemCount);
    /* Select k items with smallest vi values */
    for (i = 0; i < k; i++)
        initMedoids[i] = vIndexes[i];
    
    free(v);
    free(vIndexes);
    free(distSums);
    
    return initMedoids;
}

void SelectionSortTwoArrays(double *array, int *array2, int arrayLength)
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
                
                temp = array2[i];
                array2[i] = array2[j];
                array2[j] = temp;
                
            }
        }
    }
}

int* initializationSpreadOut(double **distanceMatrix, int itemCount, int k)
{
    int *initMedoids;
    int i, j, l, tempSum;
    double *distSquares, distSqurSum, probTemp;
    
    initMedoids = malloc(sizeof(int)*k);
    if (initMedoids == NULL)
    {
        printf("ClusterFunctions.c: initializationSpreadOut: Failed to allocate memory for initMedoids\n");
        exit(-1);
    }
    distSquares = malloc(sizeof(double)*itemCount);
    if (distSquares == NULL)
    {
        printf("ClusterFunctions.c: initializationSpreadOut: Failed to allocate memory for distSquares\n");
        exit(-1);
    }
    
    /* Choose a random medoid */
    initMedoids[0] = rand()%itemCount;
    
    for (i = 1; i < k; i++) /* We select a new medoid after each iteration */
    {
        distSqurSum = 0;
        for (j = 0; j < itemCount; j++) /* For each item */
        {
            /* Calculate squares of distances */
            distSquares[j] = distanceMatrix[initMedoids[0]][j] * distanceMatrix[initMedoids[0]][j] ;
            for (l = 1; l < i; l++) /* We're looking for the least distance among medoids */
                if (distSquares[j] > distanceMatrix[initMedoids[l]][j]*distanceMatrix[initMedoids[l]][j])
                    distSquares[j] = distanceMatrix[initMedoids[l]][j]*distanceMatrix[initMedoids[l]][j];
        
            distSqurSum += distSquares[j]; /* Least distances sum */
        }
        
        probTemp = (rand() / (double)RAND_MAX ) * distSqurSum;      /* Choose a random number as "probability" */
        tempSum = 0;
        for (j = 0; j < itemCount; j++)
        {
            tempSum += distSquares[j];
            if (tempSum >= probTemp)        /* If an item's squared distance is >= "probability" it becomes the new medoid */
                break;
        }
        
        initMedoids[i] = j;
    }
    
   // free(distSquares);
    
    return initMedoids;
}

int* assignmentPAM(double **distanceMatrix, int *medoids, int itemCount, int k, double *totalCost)
{
    int i, j;
    int *clustersArray;
    double minDistance;
    
    clustersArray = malloc(sizeof(int)*itemCount);
    if (clustersArray == NULL)
    {
        printf("ClusterFunctions.c: assignmentPAM: Failed to allocate memory for clustersArray\n");
        exit(-1);
    }
    
    *totalCost = 0;
    for (i = 0; i < itemCount; i++) /* For each item */
    {
        /* Choose a large number as min distance */
        minDistance = (double)RAND_MAX;
        for (j = 0; j < k; j++) /* For each medoid */
        {
            /* Find medoid with min distance from item */
            if (minDistance > distanceMatrix[i][medoids[j]])
            {
                minDistance = distanceMatrix[i][medoids[j]];
                clustersArray[i] = j;   /* Item will go to medoids j cluster */
            }
        }
        *totalCost += minDistance; /* Calculate total cost */
    }
    return clustersArray;
}

int* assignmentLSH(double **distanceMatrix, int **bucketIndex, int *medoids, int itemCount, double *totalCost, int k, int hashK, int L, double radius, double threshold)
{
    int i, j, m, n, *clustersArray, searchBucket, unasgCluster;
    double unasgMinDist, tempDistance;

    clustersArray = malloc(sizeof(int)*itemCount);
    if (clustersArray == NULL)
    {
        printf("ClusterFunctions.c: assignmentLSH: Failed to allocate memory for clustersArray\n");
        exit(-1);
    }
    /* Initialize clustersArray to -1 to indicate that all points are unassigned */
    for (i = 0; i < itemCount; i++)
        clustersArray[i] = -1;
    
    /* While radius < threshold = 5 * (max distance between medoids) */
    for (; radius < threshold; radius *= 2)
    {
        for (j = 0; j < L; j++) /* For each hash table */
        {
            for (m = 0; m < k; m++) /* For each medoid */
            {
                /* Find which bucket medoid is in */
                searchBucket = bucketIndex[j][medoids[m]];
                for (n = 0; n < itemCount; n++) /* For all items */
                {
                    if (bucketIndex[j][n] == searchBucket)  /* If item in same bucket as medoid */
                    {
                        if (clustersArray[n] == -1) /* If item is unassigned */
                        {
                            if (distanceMatrix[n][medoids[m]] <= radius) /* If item in radius of medoid */
                            {
                                clustersArray[n] = m;       /* Assign item to cluster */
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* For unassigned items */
    for (i = 0; i < itemCount; i++)
    {
        unasgMinDist = (double)RAND_MAX;    /* Set min distance as large number */
        if (clustersArray[i] == -1)         /* If item not assigned to cluster */
        {
            for (j = 0; j < k; j++)     /* For each medoid */
            {
                /* Find medoid with least distance from item */
                tempDistance = distanceMatrix[i][medoids[j]];
                if (tempDistance < unasgMinDist)
                {
                    unasgMinDist = tempDistance;
                    unasgCluster = medoids[j];
                }
            }
            /* Assign item to cluster with least distance from its medoid */
            clustersArray[i] = unasgCluster;
        }
    }
    /* Calculate total cost */
    *totalCost = 0.0;
    for (i = 0; i < itemCount; i++)
        *totalCost += distanceMatrix[i][clustersArray[i]];
    
    return clustersArray;
}

int* updateALaLloyds(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int assignment, int **bucketIndex, int hashK, int L, double radius, double threshold)
{
    int i, j, m, prevMedoid, *oldMedoids, isMedoid, same, iterCounter, *nonMedoids, temp;
    double oldCost, newCost;
    
    oldMedoids = malloc(sizeof(int)*k);
    if (oldMedoids == NULL)
    {
        printf("ClusterFunctions.c: updateALaLloyds: Failed to allocate memory for oldMedoids\n");
        exit(-1);
    }
    nonMedoids = malloc(sizeof(int)*(itemCount - k));
    if (nonMedoids == NULL)
    {
        printf("ClusterFunctions.c: updateALaLloyds: Failed to allocate memory for nonMedoids\n");
        exit(-1);
    }
    
    oldCost = *totalCost;
    /* Creating an array that contains all non medoid items */
    for (i = 0, j = 0; i < itemCount; i++)
    {
        isMedoid = 0;
        for (m = 0; m < k; m++)
        {
            /* Checking if item is a medoid */
            if (medoids[m] == i)
            {
                isMedoid = 1;
                break;
            }
        }
        /* If not, then put it in nonMedoids array */
        if(isMedoid==0)
        {
            nonMedoids[j]=i;
            j++;
        }
    }
    
    iterCounter = 0;
    while (iterCounter < 1000)
    {
        iterCounter++;
        for (i = 0; i < k; i++) /* Copying current medoids in oldMedoids array in case we need to revert back */
            oldMedoids[i] = medoids[i];
        
        for (i = 0; i < k; i++) /* For every medoid */
        {
            prevMedoid = medoids[i];    /* Keep old medoid */
            for (j = 0; j < itemCount-k; j++)
            {
                /* If item not in same cluster, continue */
                if (clustersArray[nonMedoids[j]] != medoids[i])
                    continue;
                /* If item in same cluster, check if it's a better medoid than the current one */
                temp = medoids[i];
                medoids[i] = nonMedoids[j];
                nonMedoids[j] = temp;
                
                /* 0 for PAM assignment, 1 for LSH assignment */
                if (assignment == 0)
                    clustersArray = assignmentPAM(distanceMatrix, medoids, itemCount, k, &newCost);
                else
                    clustersArray = assignmentLSH(distanceMatrix, bucketIndex, medoids, itemCount, &newCost, k, hashK, L, radius, threshold);
                
                if (newCost >= oldCost)
                {
                    /* If new cost is greater, revert back */
                    temp = medoids[i];
                    medoids[i]=nonMedoids[j];
                    nonMedoids[j]=temp;
                }
                else      /* Keep new medoid and update cost */
                    oldCost = newCost;
            }
        }
        /* Continue updating until oldMedoids and (new) medoids are the same */
        same = sameMedoids(oldMedoids, medoids, k);
        if (same  == 1)
            break;
    }
    *totalCost = oldCost;
    
    return clustersArray;
}

int* updateClarans(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int numLocal, int maxNeighbor, int assignment, int **bucketIndex, int hashK, int L, double radius, double threshold)
{
    int i, j, m, *currentBestMedoids, *nonMedoids, temp, *assignm, whichChanges, whichGets, isMedoid;
    double optimalCost, localCost, currentCost;
    
    currentBestMedoids = malloc(sizeof(int)*k);
    if (currentBestMedoids == NULL)
    {
        printf("ClusterFunctions.c: updateClarans: Failed to allocate memory for currentBestMedoids\n");
        exit(-1);
    }
    nonMedoids = malloc(sizeof(int)*(itemCount - k));
    if (nonMedoids == NULL)
    {
        printf("ClusterFunctions.c: updateClarans: Failed to allocate memory for nonMedoids\n");
        exit(-1);
    }
    
    optimalCost = *totalCost;
    /* Creating an array that contains all non medoid items */
    for (i = 0, j = 0; i < itemCount; i++)
    {
        isMedoid = 0;
        for (m = 0; m < k; m++)
        {
            if (medoids[m] == i)
            {
                isMedoid = 1;
                break;
            }
        }
        
        if(isMedoid==0)
        {
            nonMedoids[j]=i;
            j++;
        }
    }
    
    for (i = 0; i < numLocal; i++)
    {
        /* Randomly select k new medoids */
        currentBestMedoids =  selectKBits(itemCount, k);
        /* Make the assignment to update (local) cost */
        if (assignment == 0)
        {
            assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &localCost);
            //free(assignm);
        }
        else
        {
            assignm = assignmentLSH(distanceMatrix, bucketIndex, currentBestMedoids, itemCount, &localCost, k, hashK, L, radius, threshold);
            //free(assignm);
        }
        
        for (j = 0; j < maxNeighbor; j++)
        {
            /* Randomly select a medoid to change */
            whichChanges = rand()%k;
            /* Randomly select a non medoid item */
            whichGets = rand()%(itemCount - k);
            /* Interchange medoid with non medoid */
            temp = currentBestMedoids[whichChanges];
            currentBestMedoids[whichChanges] = nonMedoids[whichGets];
            nonMedoids[whichGets] = temp;
            /* Make the assignment to calculate new cost */
            if (assignment == 0)
            {
                assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &currentCost);
                //free(assignm);
            }
            else
            {
                assignm = assignmentLSH(distanceMatrix, bucketIndex, currentBestMedoids, itemCount, &currentCost, k, hashK, L, radius, threshold);
               // free(assignm);
            }
            
            if (currentCost < localCost)
            {
                /* If newis lower than previous cost, restart */
                j = -1;
                localCost = currentCost;
            }
            else
            {
                /* Otherwise, revert back */
                temp = currentBestMedoids[whichChanges];
                currentBestMedoids[whichChanges] = nonMedoids[whichGets];
                nonMedoids[whichGets] = temp;
            }
        }

        /* If optimal cost is greater than local cost */
        if (optimalCost > localCost)
        {
            /* Update optimal cost */
            optimalCost = localCost;
            for (j = 0; j < k; j++)     /* Update medoids to new configuration */
            {
                medoids[j] = currentBestMedoids[j];
            }
            
            if (assignment == 0)
                assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &localCost);
            else
                assignm = assignmentLSH(distanceMatrix, bucketIndex, medoids, itemCount, &localCost, k, hashK, L, radius, threshold);
        }
    }
    
    if (assignment == 0)
        clustersArray = assignmentPAM(distanceMatrix, medoids, itemCount, k, totalCost);
    else
        clustersArray = assignmentLSH(distanceMatrix, bucketIndex, medoids, itemCount, totalCost, k, hashK, L, radius, threshold);
    
    return clustersArray;
}

void medoidDistance(double **distanceMatrix, double *minDistance, double *maxDistance, int *medoids, int k)
{
    int i, j;
    double tempDistance;
    
    *minDistance = (double)RAND_MAX;    /* Set min distance to a large number */
    *maxDistance = 0.0;
    
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < k; j++)
        {
            /* Do not calculate distance from itself */
            if (i == j)
                continue;
            
            tempDistance = distanceMatrix[medoids[i]][medoids[j]];
            if (tempDistance > *maxDistance)
                *maxDistance = tempDistance;
            else if (tempDistance < *minDistance)
                *minDistance = tempDistance;
            else
                continue;
        }
    }
}

double* silhouette(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k)
{
    int i, j, m, searchCluster, counter;
    double *a, *b, *s, sum, medoidDistance;
    
    a = malloc(sizeof(double)*itemCount);
    if (a == NULL)
    {
        printf("ClusterFunctions.c: silhouette: Failed to allocate memory for a (array)\n");
        exit(-1);
    }
    b = malloc(sizeof(double)*itemCount);
    if (b == NULL)
    {
        printf("ClusterFunctions.c: silhouette: Failed to allocate memory for b (array)\n");
        exit(-1);
    }
    s = malloc(sizeof(double)*(k+1));
    if (s == NULL)
    {
        printf("ClusterFunctions.c: silhouette: Failed to allocate memory for s (array)\n");
        exit(-1);
    }
    for (i = 0; i < (k+1); i++)
        s[i] = 0.0;
    for (i = 0; i < itemCount; i++) /* For each item */
    {
        counter = 0;
        sum = 0.0;
        /* For each item, calculate sum of distances from other items in its cluster */
        for (j = 0; j <itemCount; j++)
        {
            if (clustersArray[i] == clustersArray[j])
            {
                counter++;
                sum += distanceMatrix[i][j];
            }
        }
        
        a[i] = sum/counter; /* Calculate mean distance from items in same cluster */
        counter = 0;
        sum = 0.0;
        medoidDistance = (double)RAND_MAX; /* Set min distance to neighbour medoids to large number */
        /* For each medoid, find medoid with min distance from item */
        for (m = 0; m < k; m++)
        {
            if (clustersArray[i] == m)  /* If medoid from same cluster, continue */
                continue;
            if ((medoidDistance > distanceMatrix[i][medoids[m]]) && (distanceMatrix[i][medoids[m]] != 0))
            {
                medoidDistance = distanceMatrix[i][medoids[m]];
                searchCluster = m;
            }
        }
        /* Calculate distance from items of searchCluster (cluster with min distance from item) */
        for (j = 0; j < itemCount; j++)
        {
            if (clustersArray[j] == searchCluster)
            {
                counter++;
                sum += distanceMatrix[i][j];
            }
        }
        b[i] = sum/counter; /* Calculate mean distance from items of neighbor cluster */
        
        if (b[i] > a[i])
            sum = (b[i] - a[i])/b[i];
        else
            sum = (b[i] - a[i])/a[i];
        s[i] = sum/counter;
        s[k] += sum;
    }
    
    s[k] /= itemCount;
    
    free(a);
    free(b);
    
    return s;
}

int* getClusterSize(int *clustersArray, int *medoids, int itemCount, int k)
{
    int *clusterSize, i, j;
    
    clusterSize = malloc(sizeof(int)*k);
    if (clusterSize == NULL)
    {
        printf("ClusterFunctions.c: getClusterSize: Failed to allocate memory for clusterSize\n");
        exit(-1);
    }
    
    for (i = 0; i < k; i++)
    {
        clusterSize[i] = 0;
        for (j = 0; j < itemCount; j++)
        {
            if (clustersArray[j] == i)
                clusterSize[i]++;
        }
    }
    return clusterSize;
}